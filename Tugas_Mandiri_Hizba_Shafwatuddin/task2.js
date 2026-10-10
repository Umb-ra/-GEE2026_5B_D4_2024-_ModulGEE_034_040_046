
/************************************************************
 KALIBRASI EKSPERIMENTAL NDBI
 5 kandidat terbangun + 5 kandidat bukan terbangun
 Pembanding ambang NDBI awal = -0.1
************************************************************/

// 1. ROI ASLI
var roi = ee.Geometry({
  geodesic: false,
  type: 'Polygon',
  coordinates: [[
    [106.53598136542733, -6.6221756519036195],
    [107.23910636542733, -6.6221756519036195],
    [107.23910636542733, -5.928728445061094],
    [106.53598136542733, -5.928728445061094],
    [106.53598136542733, -6.6221756519036195]
  ]]
});

// 2. LANDSAT 8 COLLECTION 2 LEVEL 2
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');

var filtered = l8
  .filterBounds(roi)
  .filterDate('2022-01-01', '2023-01-01');

print('Jumlah citra sebelum masking:', filtered.size());

// 3. MASKING AWAN DAN SKALA REFLEKTANSI
function cloudMask(image) {
  var qa = image.select('QA_PIXEL');

  var dilated = 1 << 1;
  var cirrus = 1 << 2;
  var cloud = 1 << 3;
  var shadow = 1 << 4;

  var mask = qa.bitwiseAnd(dilated).eq(0)
    .and(qa.bitwiseAnd(cirrus).eq(0))
    .and(qa.bitwiseAnd(cloud).eq(0))
    .and(qa.bitwiseAnd(shadow).eq(0));

  var scaled = image.select(
    ['SR_B2', 'SR_B3', 'SR_B4',
     'SR_B5', 'SR_B6', 'SR_B7'],
    ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']
  ).multiply(0.0000275).add(-0.2);

  return scaled.updateMask(mask)
    .copyProperties(image, image.propertyNames());
}

// 4. KOMPOSIT TAHUN 2022
var median = filtered
  .map(cloudMask)
  .median()
  .clip(roi);

// 5. HITUNG NDBI DAN NDWI
var ndbi = median.expression(
  '(SWIR - NIR) / (SWIR + NIR)', {
    SWIR: median.select('B6'),
    NIR: median.select('B5')
  }).rename('NDBI');

var ndwi = median.expression(
  '(GREEN - NIR) / (GREEN + NIR)', {
    GREEN: median.select('B3'),
    NIR: median.select('B5')
  }).rename('NDWI');

// 6. CITRA UNTUK SAMPLING
var sampleImage = ndbi.addBands(ndwi);

// Buang piksel air berdasarkan kriteria NDWI <= 0.1
var nonWater = sampleImage.updateMask(ndwi.lte(0.1));

// Ambil sampel acak untuk mencari kandidat
var samples = nonWater.sample({
  region: roi,
  scale: 30,
  numPixels: 5000,
  seed: 42,
  geometries: true,
  dropNulls: true,
  tileScale: 4
});

print('Jumlah sampel valid:', samples.size());

// 7. PILIH 5 KANDIDAT DENGAN NDBI TERTINGGI
// Ini kandidat terbangun, perlu diperiksa secara visual.
var builtCandidates = samples
  .sort('NDBI', false)
  .limit(5)
  .map(function(f) {
    return f.set('kelas_kandidat', 'Terbangun');
  });

// 8. PILIH 5 KANDIDAT DENGAN NDBI TERENDAH
var nonBuiltCandidates = samples
  .sort('NDBI', true)
  .limit(5)
  .map(function(f) {
    return f.set('kelas_kandidat', 'Bukan terbangun');
  });

// 9. TAMPILKAN KOORDINAT DAN NILAI
print('5 kandidat terbangun:', builtCandidates);
print('5 kandidat bukan terbangun:', nonBuiltCandidates);

Map.centerObject(roi, 10);

Map.addLayer(median, {
  bands: ['B4', 'B3', 'B2'],
  min: 0.02,
  max: 0.30,
  gamma: 1.2
}, 'Landsat 8 RGB');

Map.addLayer(builtCandidates, {
  color: 'ff0000'
}, 'Kandidat terbangun');

Map.addLayer(nonBuiltCandidates, {
  color: '00ffff'
}, 'Kandidat bukan terbangun');

// 10. HITUNG RATA-RATA KEDUA KELOMPOK
var meanBuilt = ee.Number(
  builtCandidates.aggregate_mean('NDBI')
);

var meanNonBuilt = ee.Number(
  nonBuiltCandidates.aggregate_mean('NDBI')
);

// Ambang baru: titik tengah kedua rata-rata
var thresholdNew = meanBuilt
  .add(meanNonBuilt)
  .divide(2);

print('Rata-rata NDBI kandidat terbangun:', meanBuilt);
print('Rata-rata NDBI kandidat bukan terbangun:', meanNonBuilt);
print('Ambang awal:', -0.1);
print('Ambang baru eksperimental:', thresholdNew);

// 11. KLASIFIKASI DENGAN AMBANG LAMA DAN BARU
var waterMask = ndwi.lte(0.1);

var builtOld = ndbi.gt(-0.1)
  .and(waterMask)
  .rename('built_old');

var builtNew = ndbi.gt(thresholdNew)
  .and(waterMask)
  .rename('built_new');

// 12. HITUNG LUAS DALAM HEKTARE
function calculateArea(binaryImage) {
  var areaImage = ee.Image.pixelArea()
    .updateMask(binaryImage)
    .rename('area');

  var result = areaImage.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 30,
    maxPixels: 1e10,
    tileScale: 4
  });

  return ee.Number(
    ee.Algorithms.If(
      result.get('area'),
      result.get('area'),
      0
    )
  ).divide(10000);
}

var areaOld = calculateArea(builtOld);
var areaNew = calculateArea(builtNew);

var difference = areaNew.subtract(areaOld);

var percentChange = ee.Algorithms.If(
  areaOld.gt(0),
  difference.divide(areaOld).multiply(100),
  null
);

print('Luas ambang -0.1 (ha):', areaOld);
print('Luas ambang baru (ha):', areaNew);
print('Perubahan luas (ha):', difference);
print('Perubahan luas (%):', percentChange);

// 13. VISUALISASI PERBANDINGAN
Map.addLayer(
  builtOld.selfMask(),
  {palette: ['ff7070']},
  'Terbangun: ambang -0.1'
);

Map.addLayer(
  builtNew.selfMask(),
  {palette: ['00ffff']},
  'Terbangun: ambang baru'
);

// Kuning = piksel yang berubah kelas
var changed = builtOld.neq(builtNew);

Map.addLayer(
  changed.selfMask(),
  {palette: ['ffff00']},
  'Perubahan klasifikasi'
);
