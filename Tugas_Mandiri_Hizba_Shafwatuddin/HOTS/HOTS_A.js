
/************************************************************
 PERBANDINGAN MASKING AWAN LANDSAT 8 — 2022
 ROI: area Jabodetabekpunjur berdasarkan koordinat script
************************************************************/

// 1. ROI
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

// 2. KOLEKSI LANDSAT 8 TAHUN 2022
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');

var filtered = l8
  .filterBounds(roi)
  .filterDate('2022-01-01', '2023-01-01');

print('Jumlah citra Landsat 8:', filtered.size());

// 3. MASKING KETAT — METODE MODUL
// Piksel ditolak jika minimal satu dari empat bit menyala.
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

  var scaled = image.select(
    ['SR_B2', 'SR_B3', 'SR_B4',
     'SR_B5', 'SR_B6', 'SR_B7'],
    ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']
  ).multiply(0.0000275).add(-0.2);

  return scaled.updateMask(mask)
    .copyProperties(image, ['system:time_start']);
}

// 4. MASKING KUSTOM
// Piksel ditolak hanya jika minimal DUA bit menyala.
function maskAwanKustom(image) {
  var qa = image.select('QA_PIXEL');

  var dilated = qa.bitwiseAnd(1 << 1).gt(0);
  var cirrus  = qa.bitwiseAnd(1 << 2).gt(0);
  var cloud   = qa.bitwiseAnd(1 << 3).gt(0);
  var shadow  = qa.bitwiseAnd(1 << 4).gt(0);

  // Konversi Boolean ke bilangan agar dapat dijumlahkan.
  var jumlahBit = dilated.toInt()
    .add(cirrus.toInt())
    .add(cloud.toInt())
    .add(shadow.toInt());

  // Pertahankan piksel jika jumlah bit aktif kurang dari 2.
  var mask = jumlahBit.lt(2);

  var scaled = image.select(
    ['SR_B2', 'SR_B3', 'SR_B4',
     'SR_B5', 'SR_B6', 'SR_B7'],
    ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']
  ).multiply(0.0000275).add(-0.2);

  return scaled.updateMask(mask)
    .copyProperties(image, ['system:time_start']);
}

// 5. KOMPOSIT MEDIAN UNTUK KEDUA METODE
var medianKetat = filtered
  .map(cloudMask)
  .median()
  .clip(roi);

var medianKustom = filtered
  .map(maskAwanKustom)
  .median()
  .clip(roi);

// 6. HITUNG JUMLAH PIXEL VALID PADA KOMPOSIT
// Menghitung lokasi piksel valid, bukan jumlah observasi
// dari seluruh citra dalam koleksi.
function hitungPikselValid(image, namaBand) {
  var hasil = image.select(namaBand).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: roi,
    scale: 30,
    maxPixels: 1e10,
    tileScale: 4
  });

  return ee.Number(hasil.get(namaBand));
}

var jumlahKetat = hitungPikselValid(medianKetat, 'B5');
var jumlahKustom = hitungPikselValid(medianKustom, 'B5');

var selisihPiksel = jumlahKustom.subtract(jumlahKetat);

var persenTambahan = ee.Algorithms.If(
  jumlahKetat.gt(0),
  selisihPiksel.divide(jumlahKetat).multiply(100),
  null
);

print('Piksel valid — masking ketat:', jumlahKetat);
print('Piksel valid — masking kustom:', jumlahKustom);
print('Tambahan piksel valid:', selisihPiksel);
print('Tambahan piksel valid (%):', persenTambahan);

// 7. VISUALISASI KOMPOSIT
var rgb = {
  bands: ['B4', 'B3', 'B2'],
  min: 0.02,
  max: 0.30,
  gamma: 1.2
};

Map.centerObject(roi, 10);

Map.addLayer(medianKetat, rgb, 'RGB — masking ketat', true);
Map.addLayer(medianKustom, rgb, 'RGB — masking kustom', false);

// 8. HITUNG NDBI DAN NDWI
function hitungIndeks(image) {
  var ndbi = image.expression(
    '(SWIR - NIR) / (SWIR + NIR)', {
      SWIR: image.select('B6'),
      NIR: image.select('B5')
    }).rename('NDBI');

  var ndwi = image.expression(
    '(GREEN - NIR) / (GREEN + NIR)', {
      GREEN: image.select('B3'),
      NIR: image.select('B5')
    }).rename('NDWI');

  return {
    ndbi: ndbi,
    ndwi: ndwi
  };
}

var indeksKetat = hitungIndeks(medianKetat);
var indeksKustom = hitungIndeks(medianKustom);

// 9. KLASIFIKASI BANGUNAN DENGAN AMBANG YANG SAMA
var bangunanKetat = indeksKetat.ndbi.gt(-0.1)
  .and(indeksKetat.ndwi.lte(0.1));

var bangunanKustom = indeksKustom.ndbi.gt(-0.1)
  .and(indeksKustom.ndwi.lte(0.1));

Map.addLayer(
  bangunanKetat.selfMask(),
  {palette: ['red']},
  'Bangunan — masking ketat',
  false
);

Map.addLayer(
  bangunanKustom.selfMask(),
  {palette: ['cyan']},
  'Bangunan — masking kustom',
  false
);

// 10. VISUALISASI PERUBAHAN MASKING
// 1 = valid pada kedua metode
// 2 = valid hanya pada metode kustom
// 3 = valid hanya pada metode ketat
var validKetat = medianKetat.select('B5').mask().gt(0);
var validKustom = medianKustom.select('B5').mask().gt(0);

var statusValid = ee.Image(0)
  .where(validKetat.and(validKustom), 1)
  .where(validKustom.and(validKetat.not()), 2)
  .where(validKetat.and(validKustom.not()), 3)
  .updateMask(validKetat.or(validKustom));

Map.addLayer(
  statusValid,
  {
    min: 1,
    max: 3,
    palette: ['808080', '00FFFF', 'FF00FF']
  },
  'Perbandingan piksel valid',
  false
);
