
/************************************************************
 EKSPERIMEN K-MEANS LANDSAT 9
 K = 5, 10, 20
 5 pengulangan setiap K
 Kandidat lahan terbangun = cluster 2 dan 3
************************************************************/

// 1. Wilayah kajian
var geometry2 = ee.Geometry.Polygon(
  [[[111.56510661311408, -7.125131151145682],
    [111.56510661311408, -7.178102498857165],
    [111.62192652888557, -7.178102498857165],
    [111.62192652888557, -7.125131151145682]]],
  null, false
);

Map.centerObject(geometry2, 12);
Map.addLayer(geometry2, {color: 'yellow'}, 'Wilayah kajian');

// 2. Ambil citra Landsat 9
var dataset = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
  .filterBounds(geometry2)
  .filterMetadata('CLOUD_COVER_LAND', 'less_than', 50)
  .sort('CLOUD_COVER_LAND');

print('Jumlah citra tersedia:', dataset.size());
print('Citra pertama:', dataset.first());

// 3. Pilih citra
var img = ee.Image(dataset.first());

// Band reflektansi
var bands = [
  'SR_B1', 'SR_B2', 'SR_B3', 'SR_B4',
  'SR_B5', 'SR_B6', 'SR_B7'
];

var img_pilih = img.select(bands).clip(geometry2);

Map.addLayer(img_pilih, {
  bands: ['SR_B4', 'SR_B3', 'SR_B2'],
  min: 7061,
  max: 14103,
  gamma: 1
}, 'Landsat 9 RGB');

// 4. Buat sampel pelatihan tetap
var data_training = img_pilih.sample({
  region: geometry2,
  scale: 30,
  numPixels: 500,
  seed: 42,
  geometries: false
});

print('Sampel pelatihan:', data_training.size());

// 5. Parameter eksperimen
var jumlahCluster = [5, 10, 20];
var jumlahUlangan = 5;

// 6. Pelatihan berulang dan perhitungan luas
var hasilEksperimen = [];

jumlahCluster.forEach(function(k) {

  for (var r = 0; r < jumlahUlangan; r++) {

    // Seed berbeda untuk inisialisasi K-means.
    // Parameter lain dipertahankan tetap.
    var seed = 100 + r;

    var clusterer = ee.Clusterer.wekaKMeans(
      k,       // jumlah cluster
      0,       // metode inisialisasi
      false,   // canopies
      100,     // maxCandidates
      10000,   // periodicPruning
      2,       // minDensity
      -1.5,    // t1
      -1.0,    // t2
      'Euclidean',
      500,     // maxIterations
      false,   // preserveOrder
      false,   // fast
      seed     // seed inisialisasi
    ).train(data_training);

    // Terapkan model pada citra
    var clustered = img_pilih.cluster(clusterer);

    // Kandidat kelas terbangun: cluster 2 atau 3
    var kandidat = clustered.eq(2)
      .or(clustered.eq(3))
      .selfMask();

    // Luas setiap piksel dalam meter persegi
    var areaImage = ee.Image.pixelArea()
      .rename('area')
      .updateMask(kandidat);

    // Jumlahkan luas kandidat dalam wilayah kajian
    var areaDict = areaImage.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: geometry2,
      scale: 30,
      maxPixels: 1e9,
      tileScale: 4
    });

    var areaM2 = ee.Number(
      ee.Algorithms.If(
        areaDict.get('area'),
        areaDict.get('area'),
        0
      )
    );

    // Simpan hasil setiap pengulangan
    hasilEksperimen.push(ee.Feature(null, {
      'K': k,
      'ulangan': r + 1,
      'seed': seed,
      'luas_m2': areaM2,
      'luas_ha': areaM2.divide(10000),
      'luas_km2': areaM2.divide(1e6)
    }));

    // Tampilkan hasil pengulangan pertama tiap K
    if (r === 0) {
      Map.addLayer(
        kandidat,
        {palette: ['ff7070']},
        'K=' + k + ', seed=' + seed
      );
    }
  }
});

// 7. Tampilkan tabel hasil
var tabel = ee.FeatureCollection(hasilEksperimen)
  .sort('K');

print('Tabel luas kandidat lahan terbangun:', tabel);

// 8. Ringkasan statistik luas untuk setiap K
var ringkasan = ee.FeatureCollection(
  jumlahCluster.map(function(k) {

    var subset = tabel.filter(ee.Filter.eq('K', k));

    var statistik = subset.reduceColumns({
      reducer: ee.Reducer.mean()
        .combine({
          reducer2: ee.Reducer.stdDev(),
          sharedInputs: true
        }),
      selectors: ['luas_ha']
    });

    return ee.Feature(null, {
      'K': k,
      'rata_luas_ha': statistik.get('mean'),
      'sd_luas_ha': statistik.get('stdDev'),
      'jumlah_ulangan': subset.size()
    });
  })
);

print('Ringkasan luas per jumlah cluster:', ringkasan);
