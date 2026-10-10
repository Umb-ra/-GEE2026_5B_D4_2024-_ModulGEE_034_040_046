
/************************************************************
 TASK III — URBAN EXPANSION KOTA BANDUNG
 Tahun: 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2024
************************************************************/

// 1. ROI: KOTA BANDUNG
// Periksa jumlah fitur dan batas yang muncul di Console/Map.
var admin = ee.FeatureCollection('FAO/GAUL/2015/level2');

var bandungFC = admin
  .filter(ee.Filter.eq('ADM0_NAME', 'Indonesia'))
  .filter(ee.Filter.eq('ADM1_NAME', 'Jawa Barat'))
  .filter(ee.Filter.eq('ADM2_NAME', 'Bandung'));

print('Jumlah fitur bernama Bandung:', bandungFC.size());
print('Nama fitur yang ditemukan:', bandungFC.aggregate_array('ADM2_NAME'));

var roi = bandungFC.geometry();

Map.centerObject(roi, 11);
Map.addLayer(
  bandungFC.style({color: 'FFFFFF', fillColor: '00000000', width: 2}),
  {},
  'Batas ROI'
);

// 2. KOLEKSI LANDSAT
var l4 = ee.ImageCollection('LANDSAT/LT04/C02/T1_L2');
var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2');
var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2');
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');

var yearList = [1990, 1995, 2000, 2005, 2010, 2015, 2020, 2024];

function filterCol(col, geometry, start, end) {
  return col
    .filterBounds(geometry)
    .filterDate(start, end);
}

// 3. CLOUD MASK LANDSAT 4/5/7
function cloudMaskTm(image) {
  var qa = image.select('QA_PIXEL');

  var mask = qa.bitwiseAnd(1 << 1).eq(0)  // Dilated cloud
    .and(qa.bitwiseAnd(1 << 3).eq(0))     // Cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0))     // Cloud shadow
    .and(qa.bitwiseAnd(1 << 5).eq(0));    // Snow

  return image.select(
    ['SR_B1', 'SR_B2', 'SR_B3',
     'SR_B4', 'SR_B5', 'SR_B7'],
    ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']
  ).updateMask(mask)
   .copyProperties(image, ['system:time_start']);
}

// 4. CLOUD MASK LANDSAT 8/9
function cloudMaskOli(image) {
  var qa = image.select('QA_PIXEL');

  var mask = qa.bitwiseAnd(1 << 1).eq(0)
    .and(qa.bitwiseAnd(1 << 2).eq(0))     // Cirrus
    .and(qa.bitwiseAnd(1 << 3).eq(0))     // Cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0))     // Cloud shadow
    .and(qa.bitwiseAnd(1 << 5).eq(0));    // Snow

  return image.select(
    ['SR_B2', 'SR_B3', 'SR_B4',
     'SR_B5', 'SR_B6', 'SR_B7'],
    ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']
  ).updateMask(mask)
   .copyProperties(image, ['system:time_start']);
}

// 5. KOMPOSIT MEDIAN
function landsat457(geometry, start, end) {
  var col = filterCol(l4, geometry, start, end)
    .merge(filterCol(l5, geometry, start, end))
    .merge(filterCol(l7, geometry, start, end));

  return col.map(cloudMaskTm).median().clip(geometry);
}

function landsat89(geometry, start, end) {
  var col = filterCol(l8, geometry, start, end)
    .merge(filterCol(l9, geometry, start, end));

  return col.map(cloudMaskOli).median().clip(geometry);
}

// 6. KLASIFIKASI DAN LUAS TAHUNAN
var builtImages = yearList.map(function(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = ee.Date.fromYMD(year + 1, 1, 1);

  var raw = (year < 2014)
    ? landsat457(roi, start, end)
    : landsat89(roi, start, end);

  // Skalakan reflektansi Collection 2
  var image = raw.multiply(0.0000275).add(-0.2);

  // RGB false color untuk inspeksi
  Map.addLayer(
    image,
    {
      bands: ['B5', 'B6', 'B2'],
      min: 0.0,
      max: 0.4,
      gamma: 1.2
    },
    'Komposit_' + year,
    false
  );

  // NDBI = (SWIR1 - NIR) / (SWIR1 + NIR)
  var ndbi = image.expression(
    '(SWIR - NIR) / (SWIR + NIR)', {
      SWIR: image.select('B6'),
      NIR: image.select('B5')
    }).rename('NDBI');

  // MNDWI = (Green - SWIR1) / (Green + SWIR1)
  var mndwi = image.expression(
    '(GREEN - SWIR) / (GREEN + SWIR)', {
      GREEN: image.select('B3'),
      SWIR: image.select('B6')
    }).rename('MNDWI');

  // Ambang eksploratif yang dipertahankan dari kode awal
  var builtMask = ndbi.gt(-0.1)
    .and(mndwi.lte(0))
    .rename('built_mask');

  // Simpan nilai tahun hanya di piksel terbangun
  var builtYear = ee.Image.constant(year)
    .toUint16()
    .rename('built')
    .updateMask(builtMask)
    .clip(roi);

  // Luas piksel terbangun dalam hektare
  var areaHa = ee.Image.pixelArea()
    .divide(10000)
    .rename('area_ha')
    .updateMask(builtMask)
    .clip(roi);

  Map.addLayer(
    builtYear,
    {min: year, max: year, palette: ['FF0000']},
    'Terbangun_' + year,
    false
  );

  return builtYear
    .addBands(areaHa)
    .set('year', year)
    .set('system:time_start', start.millis());
});

var builtCol = ee.ImageCollection.fromImages(builtImages);

// 7. PETA TAHUN PERTAMA TERDETEKSINYA LAHAN TERBANGUN
// min() mengambil tahun paling awal dari piksel yang terdeteksi
// sebagai terbangun pada salah satu tahun pengamatan.
var urbanExpansion = builtCol.select('built').min();

var palette = [
  '800080', // 1990
  '0000FF', // 1995
  '00FFFF', // 2000
  '008000', // 2005
  'FFFF00', // 2010
  'FFA500', // 2015
  'FF0000', // 2020
  'FFFFFF'  // 2024
];

// Ubah tahun menjadi kelas indeks 0–7 untuk legenda kategoris
var classValues = [0, 1, 2, 3, 4, 5, 6, 7];

var expansionClass = urbanExpansion
  .remap(yearList, classValues)
  .rename('first_detection_class');

Map.addLayer(
  expansionClass,
  {
    min: 0,
    max: 7,
    palette: palette
  },
  'Urban expansion — tahun deteksi pertama',
  true
);

// 8. HITUNG LUAS TERBANGUN PER TAHUN
var areaFeatures = yearList.map(function(year) {
  var annualImage = ee.Image(
    builtCol.filter(ee.Filter.eq('year', year)).first()
  );

  var result = annualImage.select('area_ha').reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 30,
    maxPixels: 1e10,
    tileScale: 4
  });

  var area = ee.Number(
    ee.Algorithms.If(
      result.contains('area_ha'),
      result.get('area_ha'),
      0
    )
  );

  return ee.Feature(null, {
    year: year,
    area_ha: area
  });
});

var areaFC = ee.FeatureCollection(areaFeatures);

print('Tabel luas terbangun tahunan (ha):', areaFC);

// 9. GRAFIK LUAS TERBANGUN
var chart = ui.Chart.feature.byFeature(
  areaFC.sort('year'),
  'year',
  ['area_ha']
)
.setChartType('LineChart')
.setOptions({
  title: 'Perubahan Luas Lahan Terbangun — Kota Bandung',
  hAxis: {
    title: 'Tahun',
    format: '####',
    gridlines: {count: 8}
  },
  vAxis: {
    title: 'Luas lahan terbangun (hektare)',
    minValue: 0
  },
  lineWidth: 2,
  pointSize: 5,
  legend: {position: 'none'}
});

print(chart);

// 10. LEGENDA PETA EKSPANSI
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 12px'
  }
});

legend.add(ui.Label({
  value: 'Tahun deteksi pertama',
  style: {fontWeight: 'bold', fontSize: '14px'}
}));

yearList.forEach(function(year, index) {
  var colorBox = ui.Label('', {
    backgroundColor: '#' + palette[index],
    padding: '8px',
    margin: '0 8px 4px 0',
    border: '1px solid black'
  });

  var label = ui.Label(String(year), {
    margin: '0 0 4px 0'
  });

  legend.add(ui.Panel(
    [colorBox, label],
    ui.Panel.Layout.flow('horizontal')
  ));
});

Map.add(legend);

