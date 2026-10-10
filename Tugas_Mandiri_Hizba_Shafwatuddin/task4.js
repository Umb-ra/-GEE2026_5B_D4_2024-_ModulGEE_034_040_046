function cloudMask(image) {
  var qa = image.select('QA_PIXEL');      // band penanda kualitas piksel
  var dilated = 1 << 1;                   // bit 1: awan melebar
  var cirrus  = 1 << 2;                   // bit 2: awan cirrus (tipis)
  var cloud   = 1 << 3;                   // bit 3: awan
  var shadow  = 1 << 4;                   // bit 4: bayangan awan
  var mask = qa.bitwiseAnd(dilated).eq(0)
    .and(qa.bitwiseAnd(cirrus).eq(0))
    .and(qa.bitwiseAnd(cloud).eq(0))
    .and(qa.bitwiseAnd(shadow).eq(0));
  var rename = image.select(
    ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
    ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']);
  var scaled = rename.multiply(0.0000275).add(-0.2);
  return scaled.updateMask(mask);
}


var roi = ee.Geometry({
  geodesic: false,
  type: 'Polygon',
  coordinates: [[
    [103.47443690212873, -1.7079691352593016],
    [103.7202559939256, -1.7079691352593016],
    [103.7202559939256, -1.469108505971908],
    [103.47443690212873, -1.469108505971908],
    [103.47443690212873, -1.7079691352593016]
  ]]
});

var filtered = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filterDate('2022-01-01', '2022-12-31');

// Fungsi cloudMask identik Task II -- salin ke skrip ini
var image = filtered.map(cloudMask).median().clip(roi);

Map.addLayer(image, {bands: ['B5', 'B6', 'B2'], min: [0.1, 0.05, 0], max: [0.4, 0.3, 0.15]}, 'Image');
Map.centerObject(roi, 11);

var bandMap = {
  NIR: image.select('B5'),
  SWIR: image.select('B6'),
  GREEN: image.select('B3'),
  RED: image.select('B4'),
  BLUE: image.select('B2')
};


// Badan air
var ndwi = image.expression('(GREEN - NIR) / (GREEN + NIR)', bandMap);
var air = ndwi.gt(0.1);
var bukanAir = air.eq(0);
Map.addLayer(air, {palette: ['black', 'blue'], min: 0, max: 1}, 'Air', false);

var srtm = ee.Image('USGS/SRTMGL1_003').clip(roi).updateMask(bukanAir);
Map.addLayer(srtm, {min: 0, max: 200, palette: ['green', 'yellow', 'red']}, 'Elevasi', false);

var harkatElevasi = srtm.where(srtm.gt(20), 1)
  .where(srtm.gt(15).and(srtm.lte(20)), 2)
  .where(srtm.gt(10).and(srtm.lte(15)), 3)
  .where(srtm.gt(5).and(srtm.lte(10)), 4)
  .where(srtm.lte(5), 5);
Map.addLayer(harkatElevasi, {min: 1, max: 5, palette: ['green', 'yellow', 'red']}, 'Harkat elevasi', false);

var tpi = srtm.subtract(srtm.focalMean(3))
  .reproject('EPSG:4326', null, 30)
  .updateMask(bukanAir)
  .clip(roi);
Map.addLayer(tpi, {min: -10, max: 10, palette: ['blue', 'yellow', 'red']}, 'TPI', false);

var harkatTPI = tpi.where(tpi.gt(0), 1)
  .where(tpi.gt(-0.25).and(tpi.lte(0)), 2)
  .where(tpi.gt(-0.5).and(tpi.lte(-0.25)), 3)
  .where(tpi.gt(-0.75).and(tpi.lte(-0.5)), 4)
  .where(tpi.lte(-0.75), 5);
Map.addLayer(harkatTPI, {min: 1, max: 5, palette: ['green', 'yellow', 'red']}, 'Harkat TPI', false);

var jarakAir = air.fastDistanceTransform()
  .sqrt().multiply(30)      // jarak kuadrat-piksel -> meter
  .updateMask(bukanAir)
  .clip(roi)
  .reproject('EPSG:4326', null, 30);
Map.addLayer(jarakAir, {min: 0, max: 1000, palette: ['red', 'yellow', 'green']}, 'Jarak dari air', false);

// Harkat jarak (Tabel 4)
var harkatJarakAir = jarakAir.where(jarakAir.gt(200), 1)
  .where(jarakAir.gt(150).and(jarakAir.lte(200)), 2)
  .where(jarakAir.gt(100).and(jarakAir.lte(150)), 3)
  .where(jarakAir.gt(50).and(jarakAir.lte(100)), 4)
  .where(jarakAir.lte(50), 5);
Map.addLayer(harkatJarakAir, {min: 1, max: 5, palette: ['green', 'yellow', 'red']}, 'Harkat jarak dari air', false);

var ndvi = image.expression('(NIR - RED) / (NIR + RED)', bandMap).updateMask(bukanAir);
Map.addLayer(ndvi, {min: -1, max: 1, palette: ['blue', 'white', 'green']}, 'NDVI', false);

var harkatNDVI = ndvi.where(ndvi.gt(0.5), 1)
  .where(ndvi.gt(0.25).and(ndvi.lte(0.5)), 2)
  .where(ndvi.gt(0).and(ndvi.lte(0.25)), 3)
  .where(ndvi.gt(-0.25).and(ndvi.lte(0)), 4)
  .where(ndvi.lte(-0.25), 5);
Map.addLayer(harkatNDVI, {min: 1, max: 5, palette: ['green', 'yellow', 'red']}, 'Harkat NDVI', false);

var kerawanan = harkatJarakAir.add(harkatElevasi).add(harkatTPI).add(harkatNDVI);
Map.addLayer(kerawanan, {min: 1, max: 20, palette: ['green', 'yellow', 'red']}, 'Kerawanan', false);

var harkatKerawanan = kerawanan.where(kerawanan.gt(16), 5)
  .where(kerawanan.gt(12).and(kerawanan.lte(16)), 4)
  .where(kerawanan.gt(8).and(kerawanan.lte(12)), 3)
  .where(kerawanan.gt(4).and(kerawanan.lte(8)), 2)
  .where(kerawanan.lte(4), 1);
Map.addLayer(harkatKerawanan, {min: 1, max: 5, palette: ['green', 'yellow', 'red']}, 'Harkat kerawanan', false);

var penduduk = ee.ImageCollection('JRC/GHSL/P2023A/GHS_POP')
  .sort('system:time_start', false)
  .first().updateMask(bukanAir).clip(roi);
Map.addLayer(penduduk, {min: 0, max: 250, palette: ['black', 'red', 'white']}, 'Penduduk', false);

var harkatPenduduk = penduduk.where(penduduk.gt(300), 5)
  .where(penduduk.gt(200).and(penduduk.lte(300)), 4)
  .where(penduduk.gt(100).and(penduduk.lte(200)), 3)
  .where(penduduk.gt(10).and(penduduk.lte(100)), 2)
  .where(penduduk.lte(10), 1);
Map.addLayer(harkatPenduduk, {min: 1, max: 5, palette: ['green', 'yellow', 'red']}, 'Harkat penduduk', false);


var wc = ee.ImageCollection('ESA/WorldCover/v200').first().clip(roi).updateMask(bukanAir);

var lahanTerbangun = wc.eq(50);   // kelas 50 = built-up
var pertanian = wc.eq(40);        // kelas 40 = cropland

var area = ee.Image.pixelArea();

// Luas terbangun per piksel 30 m (m persegi)
var lahanTerbangunArea = lahanTerbangun.multiply(area)
  .reduceResolution(ee.Reducer.sum(), true)
  .reproject('EPSG:4326', null, 30);
Map.addLayer(lahanTerbangunArea, {min: 0, max: 100, palette: ['black', 'red', 'white']}, 'Luas lahan terbangun', false);

var harkatLahanTerbangun = lahanTerbangunArea.where(lahanTerbangunArea.gt(100), 5)
  .where(lahanTerbangunArea.gt(80).and(lahanTerbangunArea.lte(100)), 4)
  .where(lahanTerbangunArea.gt(60).and(lahanTerbangunArea.lte(80)), 3)
  .where(lahanTerbangunArea.gt(40).and(lahanTerbangunArea.lte(60)), 2)
  .where(lahanTerbangunArea.lte(40), 1);
Map.addLayer(harkatLahanTerbangun, {min: 1, max: 5, palette: ['green', 'yellow', 'red']}, 'Harkat lahan terbangun', false);

// Luas pertanian per piksel 30 m + harkat
var pertanianArea = pertanian.multiply(area)
  .reduceResolution(ee.Reducer.sum(), true)
  .reproject('EPSG:4326', null, 30);
Map.addLayer(pertanianArea, {min: 0, max: 100, palette: ['black', 'lightgreen', 'darkgreen']}, 'Luas pertanian', false);

var harkatPertanian = pertanianArea.where(pertanianArea.gt(100), 5)
  .where(pertanianArea.gt(80).and(pertanianArea.lte(100)), 4)
  .where(pertanianArea.gt(60).and(pertanianArea.lte(80)), 3)
  .where(pertanianArea.gt(40).and(pertanianArea.lte(60)), 2)
  .where(pertanianArea.lte(40), 1);
Map.addLayer(harkatPertanian, {min: 1, max: 5, palette: ['green', 'yellow', 'red']}, 'Harkat pertanian', false);


var eksposur = harkatPenduduk.add(harkatLahanTerbangun).add(harkatPertanian);
Map.addLayer(eksposur, {min: 1, max: 15, palette: ['green', 'yellow', 'red']}, 'Eksposur', false);

var harkatEksposur = eksposur.where(eksposur.gt(12), 5)
  .where(eksposur.gt(9).and(eksposur.lte(12)), 4)
  .where(eksposur.gt(6).and(eksposur.lte(9)), 3)
  .where(eksposur.gt(3).and(eksposur.lte(6)), 2)
  .where(eksposur.lte(3), 1);
Map.addLayer(harkatEksposur, {min: 1, max: 5, palette: ['green', 'yellow', 'red']}, 'Harkat eksposur', false);


// PDB per piksel (USD)
var pdb = ee.Image('projects/sat-io/open-datasets/GRIDDED_HDI_GDP/GDP_PPP_1990_2015_5arcmin_v2')
  .select('b26').clip(roi).updateMask(bukanAir);
Map.addLayer(pdb, {min: 0, max: 1e8, palette: ['black', 'green', 'white']}, 'PDB', false);

var harkatPDB = pdb.where(pdb.gt(40 * 1e6), 1)
  .where(pdb.gt(30 * 1e6).and(pdb.lte(40 * 1e6)), 2)
  .where(pdb.gt(20 * 1e6).and(pdb.lte(30 * 1e6)), 3)
  .where(pdb.gt(10 * 1e6).and(pdb.lte(20 * 1e6)), 4)
  .where(pdb.lte(10 * 1e6), 5);
Map.addLayer(harkatPDB, {min: 1, max: 5, palette: ['green', 'yellow', 'red']}, 'Harkat PDB', false);

// IPM (0--1)
var ipm = ee.Image('projects/sat-io/open-datasets/GRIDDED_HDI_GDP/HDI_1990_2015_v2')
  .select('b26').clip(roi).updateMask(bukanAir);
Map.addLayer(ipm, {min: 0.5, max: 1, palette: ['black', 'blue', 'white']}, 'IPM', false);

var harkatIPM = ipm.where(ipm.gt(0.85), 1)
  .where(ipm.gt(0.8).and(ipm.lte(0.85)), 2)
  .where(ipm.gt(0.75).and(ipm.lte(0.8)), 3)
  .where(ipm.gt(0.7).and(ipm.lte(0.75)), 4)
  .where(ipm.lte(0.7), 5);
Map.addLayer(harkatIPM, {min: 1, max: 5, palette: ['green', 'yellow', 'red']}, 'Harkat IPM', false);

var kerentanan = harkatPDB.add(harkatIPM);
Map.addLayer(kerentanan, {min: 1, max: 10, palette: ['green', 'yellow', 'red']}, 'Kerentanan', false);

var harkatKerentanan = kerentanan.where(kerentanan.gt(8), 5)
  .where(kerentanan.gt(6).and(kerentanan.lte(8)), 4)
  .where(kerentanan.gt(4).and(kerentanan.lte(6)), 3)
  .where(kerentanan.gt(2).and(kerentanan.lte(4)), 2)
  .where(kerentanan.lte(2), 1);
Map.addLayer(harkatKerentanan, {min: 1, max: 5, palette: ['green', 'yellow', 'red']}, 'Harkat kerentanan', false);


var risiko = harkatKerawanan.add(harkatEksposur);   // .add(harkatKerentanan);
Map.addLayer(risiko, {min: 1, max: 10, palette: ['green', 'yellow', 'red']}, 'Risiko banjir', false);

var harkatRisiko = risiko.where(risiko.gt(8), 5)
  .where(risiko.gt(6).and(risiko.lte(8)), 4)
  .where(risiko.gt(4).and(risiko.lte(6)), 3)
  .where(risiko.gt(2).and(risiko.lte(4)), 2)
  .where(risiko.lte(2), 1);
Map.addLayer(harkatRisiko, {min: 1, max: 5, palette: ['green', 'yellow', 'red']}, 'Harkat risiko', false);


/************************************************************
 * PERBANDINGAN RISIKO BANJIR SEBELUM DAN SESUDAH KERENTANAN
 ************************************************************/

// ==========================================================
// A. RISIKO SEBELUM: KERAWANAN + EKSPOSUR
// ==========================================================

var risikoSebelum = harkatKerawanan
  .add(harkatEksposur)
  .rename('risiko_sebelum');

var harkatRisikoSebelum = risikoSebelum
  .where(risikoSebelum.gt(8), 5)
  .where(risikoSebelum.gt(6).and(risikoSebelum.lte(8)), 4)
  .where(risikoSebelum.gt(4).and(risikoSebelum.lte(6)), 3)
  .where(risikoSebelum.gt(2).and(risikoSebelum.lte(4)), 2)
  .where(risikoSebelum.lte(2), 1)
  .rename('kelas_sebelum');

Map.addLayer(
  harkatRisikoSebelum,
  {min: 1, max: 5,
   palette: ['006400', '7FFF00', 'FFFF00', 'FFA500', 'FF0000']},
  'Harkat risiko sebelum kerentanan',
  true
);


// ==========================================================
// B. RISIKO SESUDAH: KERAWANAN + EKSPOSUR + KERENTANAN
// ==========================================================

var risikoSesudah = harkatKerawanan
  .add(harkatEksposur)
  .add(harkatKerentanan)
  .rename('risiko_sesudah');

var harkatRisikoSesudah = risikoSesudah
  .where(risikoSesudah.gt(8), 5)
  .where(risikoSesudah.gt(6).and(risikoSesudah.lte(8)), 4)
  .where(risikoSesudah.gt(4).and(risikoSesudah.lte(6)), 3)
  .where(risikoSesudah.gt(2).and(risikoSesudah.lte(4)), 2)
  .where(risikoSesudah.lte(2), 1)
  .rename('kelas_sesudah');

Map.addLayer(
  harkatRisikoSesudah,
  {min: 1, max: 5,
   palette: ['006400', '7FFF00', 'FFFF00', 'FFA500', 'FF0000']},
  'Harkat risiko sesudah kerentanan',
  false
);


// ==========================================================
// C. PETA PERUBAHAN KELAS RISIKO
// ==========================================================

// Nilai positif: kelas meningkat
// Nilai nol: kelas tidak berubah
// Nilai negatif: kelas menurun

var perubahanKelas = harkatRisikoSesudah
  .subtract(harkatRisikoSebelum)
  .rename('perubahan_kelas');

Map.addLayer(
  perubahanKelas,
  {min: -4, max: 4,
   palette: ['0000FF', 'FFFFFF', 'FF0000']},
  'Perubahan kelas risiko',
  false
);


// ==========================================================
// D. FUNGSI MENGHITUNG LUAS KELAS RISIKO (HEKTARE)
// Tanpa Reducer.group dan List.filter
// ==========================================================

var luasPikselHa = ee.Image.pixelArea()
  .divide(10000)
  .rename('luas_ha');

function hitungLuasKelas(image, nama) {

  var tabel = ee.FeatureCollection(
    ee.List.sequence(1, 5).map(function(kelas) {

      // Ambil hanya piksel yang memiliki kelas tersebut
      var maskKelas = image.eq(ee.Number(kelas));

      var hasil = luasPikselHa
        .updateMask(maskKelas)
        .reduceRegion({
          reducer: ee.Reducer.sum(),
          geometry: roi,
          scale: 30,
          crs: 'EPSG:32748',
          maxPixels: 1e13,
          tileScale: 4
        });

      var luas = ee.Number(
        ee.Algorithms.If(
          hasil.contains('luas_ha'),
          hasil.get('luas_ha'),
          0
        )
      );

      return ee.Feature(null, {
        skenario: nama,
        kelas_risiko: kelas,
        luas_ha: luas
      });
    })
  );

  return tabel;
}


// ==========================================================
// E. HITUNG LUAS SEBELUM DAN SESUDAH
// ==========================================================

var luasSebelum = hitungLuasKelas(
  harkatRisikoSebelum,
  'Sebelum'
);

var luasSesudah = hitungLuasKelas(
  harkatRisikoSesudah,
  'Sesudah'
);

print(
  'Luas tiap kelas sebelum (ha)',
  luasSebelum
);

print(
  'Luas tiap kelas sesudah (ha)',
  luasSesudah
);


// ==========================================================
// F. TABEL PERBANDINGAN LUAS
// ==========================================================

var perbandinganLuas = ee.FeatureCollection(
  ee.List.sequence(1, 5).map(function(kelas) {

    var sebelum = ee.Number(
      luasSebelum
        .filter(ee.Filter.eq('kelas_risiko', kelas))
        .first()
        .get('luas_ha')
    );

    var sesudah = ee.Number(
      luasSesudah
        .filter(ee.Filter.eq('kelas_risiko', kelas))
        .first()
        .get('luas_ha')
    );

    var perubahan = sesudah.subtract(sebelum);

    var persen = ee.Algorithms.If(
      sebelum.gt(0),
      perubahan.divide(sebelum).multiply(100),
      null
    );

    return ee.Feature(null, {
      kelas_risiko: kelas,
      sebelum_ha: sebelum,
      sesudah_ha: sesudah,
      perubahan_ha: perubahan,
      perubahan_persen: persen
    });
  })
);

print(
  'Perbandingan luas kelas risiko',
  perbandinganLuas
);


// ==========================================================
// G. MATRIKS PERPINDAHAN KELAS RISIKO (HEKTARE)
// Baris = kelas sebelum
// Kolom = kelas sesudah
// ==========================================================

var matriksPerubahan = ee.FeatureCollection(
  ee.List.sequence(1, 5).map(function(kelasAwal) {

    var fiturBaris = ee.FeatureCollection(
      ee.List.sequence(1, 5).map(function(kelasAkhir) {

        var maskTransisi = harkatRisikoSebelum
          .eq(ee.Number(kelasAwal))
          .and(
            harkatRisikoSesudah.eq(ee.Number(kelasAkhir))
          );

        var hasil = luasPikselHa
          .updateMask(maskTransisi)
          .reduceRegion({
            reducer: ee.Reducer.sum(),
            geometry: roi,
            scale: 30,
            crs: 'EPSG:32748',
            maxPixels: 1e13,
            tileScale: 4
          });

        var luas = ee.Number(
          ee.Algorithms.If(
            hasil.contains('luas_ha'),
            hasil.get('luas_ha'),
            0
          )
        );

        return ee.Feature(null, {
          kelas_sebelum: kelasAwal,
          kelas_sesudah: kelasAkhir,
          luas_ha: luas
        });
      })
    );

    return fiturBaris;
  })
).flatten();

print(
  'Matriks perpindahan kelas risiko (ha)',
  matriksPerubahan
);



