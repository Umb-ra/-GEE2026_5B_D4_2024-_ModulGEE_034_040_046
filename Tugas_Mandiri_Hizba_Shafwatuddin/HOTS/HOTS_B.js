
// =====================================================
// PERBANDINGAN NDWI DAN MNDWI - LANDSAT 8 TAHUN 2022
// ROI: wilayah Jakarta dan sekitarnya
// Ambang air: indeks > 0.1
// Luas dihitung menggunakan ee.Reducer.sum()
// =====================================================

// 1. AREA PENELITIAN
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

Map.centerObject(roi, 10);
Map.addLayer(roi, {color: 'white'}, 'Batas ROI');

// 2. KOLEKSI LANDSAT 8
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');

var filtered = l8
  .filterBounds(roi)
  .filterDate('2022-01-01', '2023-01-01');

// 3. MASKING AWAN DAN SKALA REFLEKTANSI
function cloudMask(image) {
  var qa = image.select('QA_PIXEL');

  var dilated = 1 << 1;
  var cirrus  = 1 << 2;
  var cloud   = 1 << 3;
  var shadow  = 1 << 4;

  var mask = qa.bitwiseAnd(dilated).eq(0)
    .and(qa.bitwiseAnd(cirrus).eq(0))
    .and(qa.bitwiseAnd(cloud).eq(0))
    .and(qa.bitwiseAnd(shadow).eq(0));

  var bands = image.select(
    ['SR_B3', 'SR_B5', 'SR_B6'],
    ['GREEN', 'NIR', 'SWIR']
  );

  var scaled = bands.multiply(0.0000275).add(-0.2);

  return scaled.updateMask(mask)
    .copyProperties(image, ['system:time_start']);
}

// Komposit median yang sama digunakan oleh kedua indeks
var median = filtered
  .map(cloudMask)
  .median()
  .clip(roi);

print('Jumlah citra Landsat 8:', filtered.size());

// 4. MENGHITUNG NDWI DAN MNDWI
// NDWI  = (Green - NIR) / (Green + NIR)
// MNDWI = (Green - SWIR) / (Green + SWIR)

var ndwi = median.expression(
  '(GREEN - NIR) / (GREEN + NIR)', {
    GREEN: median.select('GREEN'),
    NIR: median.select('NIR')
  }
).rename('NDWI');

var mndwi = median.expression(
  '(GREEN - SWIR) / (GREEN + SWIR)', {
    GREEN: median.select('GREEN'),
    SWIR: median.select('SWIR')
  }
).rename('MNDWI');

// 5. KLASIFIKASI AIR DENGAN AMBANG YANG SAMA
var threshold = 0.1;

var waterNDWI = ndwi.gt(threshold).rename('Air_NDWI');
var waterMNDWI = mndwi.gt(threshold).rename('Air_MNDWI');

// 6. MENGHITUNG LUAS AIR (HEKTARE)
// Pixel area dikonversi dari m² ke hektare.
// Reducer.sum menjumlahkan luas semua piksel yang
// memenuhi kriteria air.

var pixelAreaHa = ee.Image.pixelArea()
  .divide(10000)
  .rename('area_ha');

var areaNDWI = pixelAreaHa
  .updateMask(waterNDWI)
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 30,
    maxPixels: 1e13,
    tileScale: 4
  });

var areaMNDWI = pixelAreaHa
  .updateMask(waterMNDWI)
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 30,
    maxPixels: 1e13,
    tileScale: 4
  });

var luasNDWI = ee.Number(areaNDWI.get('area_ha'));
var luasMNDWI = ee.Number(areaMNDWI.get('area_ha'));

print('Luas air NDWI (Ha):', luasNDWI);
print('Luas air MNDWI (Ha):', luasMNDWI);
print('Selisih luas bertanda (MNDWI - NDWI), Ha:',
      luasMNDWI.subtract(luasNDWI));
print('Selisih luas absolut (Ha):',
      luasMNDWI.subtract(luasNDWI).abs());

// 7. MEMBANDINGKAN KLASIFIKASI PER PIKSEL
// Nilai:
//  1 = hanya NDWI mengidentifikasi air
// -1 = hanya MNDWI mengidentifikasi air
//  0 = kedua indeks sepakat
//
// Perbandingan hanya dilakukan pada piksel yang valid
// pada kedua indeks.

var difference = waterNDWI
  .subtract(waterMNDWI)
  .rename('Perbedaan');

var disagreement = difference.neq(0)
  .selfMask()
  .rename('Ketidaksepakatan');

Map.addLayer(ndwi, {
  min: -1, max: 1,
  palette: ['brown', 'white', 'blue']
}, 'Nilai NDWI', false);

Map.addLayer(mndwi, {
  min: -1, max: 1,
  palette: ['brown', 'white', 'blue']
}, 'Nilai MNDWI', false);

Map.addLayer(waterNDWI.selfMask(), {
  palette: ['0000FF']
}, 'Air menurut NDWI');

Map.addLayer(waterMNDWI.selfMask(), {
  palette: ['00FFFF']
}, 'Air menurut MNDWI', false);

// Merah = hanya NDWI mendeteksi air
// Cyan = hanya MNDWI mendeteksi air
Map.addLayer(
  difference.updateMask(difference.neq(0)),
  {
    min: -1,
    max: 1,
    palette: ['00FFFF', 'FF0000']
  },
  'Peta selisih NDWI-MNDWI'
);

// 8. IDENTIFIKASI TUTUPAN LAHAN DI WILAYAH BERBEDA
// ESA WorldCover 2021
var worldcover = ee.ImageCollection(
  'ESA/WorldCover/v200'
).first().select('Map').clip(roi);

// Luas wilayah ketidaksepakatan per kelas tutupan lahan.
// Band pertama = luas (Ha), band kedua = kode kelas.
// Grouped reducer menjumlahkan luas untuk setiap kelas.

var areaByLandcover = pixelAreaHa
  .updateMask(difference.neq(0))
  .addBands(worldcover.rename('kelas'));

var landcoverStats = areaByLandcover.reduceRegion({
  reducer: ee.Reducer.sum().group({
    groupField: 1,
    groupName: 'kode_tutupan_lahan'
  }),
  geometry: roi,
  scale: 30,
  maxPixels: 1e13,
  tileScale: 4
});

print(
  'Luas ketidaksepakatan per kelas WorldCover (Ha):',
  landcoverStats
);

// Kode kelas ESA WorldCover:
// 10 = pohon
// 20 = semak
// 30 = vegetasi herba
// 40 = lahan pertanian
// 50 = kawasan terbangun
// 60 = lahan terbuka / jarang bervegetasi
// 70 = salju/es
// 80 = badan air permanen
// 90 = vegetasi akuatik
// 95 = mangrove
// 100 = lumut dan liken

// Tampilkan tutupan lahan hanya di lokasi berbeda
Map.addLayer(
  worldcover.updateMask(difference.neq(0)),
  {
    min: 10,
    max: 100,
    palette: [
      '006400', // 10 pohon
      'ffbb22', // 20 semak
      'ffff4c', // 30 vegetasi herba
      'f096ff', // 40 pertanian
      'fa0000', // 50 terbangun
      'b4b4b4', // 60 lahan terbuka
      'f0f0f0', // 70 salju/es
      '0064c8', // 80 badan air
      '0096a0', // 90 vegetasi akuatik
      '00cf75', // 95 mangrove
      'fae6a0'  // 100 lumut/liken
    ]
  },
  'Tutupan lahan pada wilayah berbeda',
  false
);

// 9. RINGKASAN HASIL DALAM TABEL
var comparison = ee.FeatureCollection([
  ee.Feature(null, {
    'Indeks': 'NDWI',
    'Ambang': threshold,
    'Luas_air_Ha': luasNDWI
  }),
  ee.Feature(null, {
    'Indeks': 'MNDWI',
    'Ambang': threshold,
    'Luas_air_Ha': luasMNDWI
  }),
  ee.Feature(null, {
    'Indeks': 'Selisih absolut',
    'Ambang': threshold,
    'Luas_air_Ha': luasNDWI.subtract(luasMNDWI).abs()
  })
]);

print('Tabel perbandingan luas air:', comparison);

