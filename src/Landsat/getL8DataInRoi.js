var place_key = "A0";
var year = '2015';
var startDateTime = year + '-06-01T10:00:00'; // Start date and time in UTC
var endDateTime = year + '-09-01T15:00:00';   // End date and time in UTC

var geom = geometry;
var bb_area_meters = geom.area({maxError: 1});
var bb_area_miles = bb_area_meters.divide(2589988.11);
var bb_pts = geom.coordinates();
print('ROI Square Miles:', bb_area_miles);

var bb_pts_coords = bb_pts.getInfo();
var roi = ee.Geometry.BBox(
  bb_pts_coords[0][3][0], bb_pts_coords[0][3][1],
  bb_pts_coords[0][1][0], bb_pts_coords[0][1][1]
);

function prepSrL8(image) {
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);

  var getFactorImg = function(factorNames) {
    var factorList = image.toDictionary().select(factorNames).values();
    return ee.Image.constant(factorList);
  };
  var scaleImg = getFactorImg([
    'REFLECTANCE_MULT_BAND_.',
    'TEMPERATURE_MULT_BAND_ST_B10']);
  var offsetImg = getFactorImg([
    'REFLECTANCE_ADD_BAND_.',
    'TEMPERATURE_ADD_BAND_ST_B10']);
    
  var scaled = image.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10'])
                    .multiply(scaleImg).add(offsetImg);
                    
  return image.addBands(scaled, null, true)
    .updateMask(qaMask).updateMask(saturationMask);
}

var l8_filtered = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filterDate(startDateTime, endDateTime)
  .map(prepSrL8)
  .select('SR.*', 'ST.*')
  .median();

var l8_clip = l8_filtered.clip(roi);

var l8_raw_filtered = ee.ImageCollection('LANDSAT/LC08/C02/T1')
  .filterDate(startDateTime, endDateTime);

var composite = ee.Algorithms.Landsat.simpleComposite({
  collection: l8_raw_filtered,
  asFloat: true
});

var l8_composite_clip = composite.clip(roi);

var thermalBands = l8_clip.select('ST_B10').toFloat(); // Convert to Float32
var opticalBands = l8_clip.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']).toFloat(); // Convert to Float32
var l8_composite_oli = l8_composite_clip.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7']).toFloat(); // Convert to Float32

var thermalBands_converted = thermalBands.expression(
  '(dn - 273.15) * 1.8 + 32.0', {
    'dn': thermalBands
  }
).rename('thermalBands_converted').toFloat(); // Ensure converted band is Float32

// Combine the images into one image
var combinedImage = ee.Image.cat([
  l8_composite_oli,
  thermalBands_converted
]);

// Print band names for verification
print('Bands in combined image:', combinedImage.bandNames());

// Display the combined image
Map.addLayer(combinedImage, {
  bands: ['B4', 'B3', 'B2'], // Example visualization for the RGB bands
  min: 0,
  max: 0.3,
  gamma: 1.4
}, 'Combined Image');

// Add the converted thermal bands to the map with a different visualization
Map.addLayer(thermalBands_converted, {palette:['purple', 'blue', 'yellow', 'orange', 'red'], min: 75, max: 155}, 'thermalBands_converted');

// Export the combined image
Export.image.toDrive({
  image: combinedImage,
  description: place_key + '_' + year,
  folder: 'L8_ARD',
  region: roi,
  scale: 30,
  fileFormat: 'GeoTIFF'
});
