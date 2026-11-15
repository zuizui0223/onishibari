// GloPL site × MODIS greenness/productivity environmental context (UNEXECUTED).
// Google Earth Engine Code Editor; requires an authenticated EE project & uploaded
// site-coordinate CSV from GloPL WorldClim pilot. No NDVI or GPP values are
// claimed until the export tasks have completed and output QC is validated.
//
// This script calculates 2001-2015 CONTEXTUAL vegetation climatology, not
// contemporaneous NDVI for original 1981-2015 pollen experiments (GloPL Year
// is PUBLICATION year and does not identify the actual experiment year).
// 250m NDVI vs 500m GPP are different native resolutions.
// Upload site_coordinates_for_ndvi_sampling.csv as point features using
// longitude and latitude columns; replace the placeholder asset ID.
var sourceSites = ee.FeatureCollection('projects/REPLACE_ME/assets/glopl_sites');

var sites = sourceSites.map(function(f) {
  var lon = ee.Number(f.get('longitude'));
  var lat = ee.Number(f.get('latitude'));
  var geom = ee.Geometry.Point([lon, lat]);
  return ee.Feature(geom, {
    site_id: f.get('site_id'),
    n_effect_size_cases: f.get('n_effect_size_cases'),
    contemporaneous_with_GloPL: false
  });
});
var start = '2001-01-01';
var endExclusive = '2016-01-01';

var indices = ee.ImageCollection('MODIS/061/MOD13Q1')
  .filterDate(start, endExclusive)
  .map(function(im) {
    var good = im.select('SummaryQA').lte(1);
    return im.select(['NDVI', 'EVI'])
      .multiply(0.0001)
      .updateMask(good)
      .copyProperties(im, ['system:time_start']);
  });

var climatologicalMonthly = ee.ImageCollection.fromImages(
  ee.List.sequence(1, 12).map(function(month) {
    return indices.filter(ee.Filter.calendarRange(month, month, 'month'))
      .mean()
      .set('month', month);
  })
);

// Captures mean greenness plus within-year cycle (not pollinator or flower counts).
var ndviMean = indices.select('NDVI').mean().rename('NDVI_mean_2001_2015');
var eviMean = indices.select('EVI').mean().rename('EVI_mean_2001_2015');
var ndviAmplitude = climatologicalMonthly.select('NDVI').max()
  .subtract(climatologicalMonthly.select('NDVI').min())
  .rename('NDVI_monthly_amplitude_2001_2015');
var eviAmplitude = climatologicalMonthly.select('EVI').max()
  .subtract(climatologicalMonthly.select('EVI').min())
  .rename('EVI_monthly_amplitude_2001_2015');

// MODIS MOD17 is 8-day GPP SUM per image, kgC/m2 after scaling 0.0001.
// Each calendar year is summed; 15 annual totals are averaged (not summed).
var gpp = ee.ImageCollection('MODIS/061/MOD17A2HGF')
  .filterDate(start, endExclusive)
  .map(function(im) {
    return im.select('Gpp').multiply(0.0001)
      .copyProperties(im, ['system:time_start']);
  });
var annualGpp = ee.ImageCollection.fromImages(
  ee.List.sequence(2001, 2015).map(function(y) {
    y = ee.Number(y);
    return gpp.filterDate(ee.Date.fromYMD(y, 1, 1),
                          ee.Date.fromYMD(y.add(1), 1, 1))
      .sum().rename('GPP_annual_kgC_m2');
  })
).mean().rename('GPP_annual_mean_kgC_m2_2001_2015');

var environmentalImage = ee.Image.cat(
  [ndviMean, eviMean, ndviAmplitude, eviAmplitude, annualGpp]
);
var selectors = ['site_id', 'n_effect_size_cases', 'contemporaneous_with_GloPL',
  'NDVI_mean_2001_2015', 'EVI_mean_2001_2015',
  'NDVI_monthly_amplitude_2001_2015',
  'EVI_monthly_amplitude_2001_2015',
  'GPP_annual_mean_kgC_m2_2001_2015'];

[250, 1000, 5000].forEach(function(radius) {
  // Simple landscape buffers; bee species-specific foraging ranges must
  // be modeled separately instead of declaring one buffer biologically correct.
  var circles = sites.map(function(f) {
    return ee.Feature(f.geometry().buffer(radius), f.toDictionary());
  });
  var sampled = environmentalImage.reduceRegions({
    collection: circles,
    reducer: ee.Reducer.mean(),
    scale: 250,
    tileScale: 8
  });
  Export.table.toDrive({
    collection: sampled,
    description: 'glopl_modis_contextual_2001_2015_' + radius + 'm',
    fileFormat: 'CSV',
    selectors: selectors
  });
});

// DO NOT merge GEDI 2019–2023 biomass with 1981–2015 pollen experiments as if
// contemporaneous. For newer events, use GEDI L4B uncertainty-aware AGB and
// an explicit year-matching / environmental stability rule.
// Direct links:
// https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13Q1
// https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD17A2HGF
// https://developers.google.com/earth-engine/datasets/catalog/LARSE_GEDI_GEDI04_B_002
