var place_key = "trans_pecos";
var year = '2018';
var startDateTime = year+'-06-01T10:00:00'; // Start date and time in UTC
var endDateTime = year+'-09-01T15:00:00';   // End date and time in UTC
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
var geom = geometry;
var bb_area_meters = geom.area({maxError: 1});
var bb_area_miles = bb_area_meters.divide(2589988.11);
var bb_pts = geom.coordinates();
print(bb_area_miles)
//////////////////////////////////////////////////////////////////////////////////////
var bb_pts_coords = bb_pts.getInfo();
print(bb_pts_coords[0][1]);

var bb_pt_a = ee.Geometry.Point([bb_pts_coords[0][3][0], bb_pts_coords[0][3][1]]);
Map.addLayer(bb_pt_a, {color: 'red'}, 'Point');

var bb_pt_b = ee.Geometry.Point([bb_pts_coords[0][1][0], bb_pts_coords[0][1][1]]);
Map.addLayer(bb_pt_b, {color: 'red'}, 'Point');

var mapCenterPoint = geom.centroid({'maxError': 1});
var mapCenterCoords = mapCenterPoint.coordinates();
var mapCenterCoords = mapCenterCoords.getInfo();
Map.setCenter(mapCenterCoords[0], mapCenterCoords[1], 8);
//////////////////////////////////////////////////////////////////////////////////////
var roi = ee.Geometry.BBox(
  bb_pts_coords[0][3][0],bb_pts_coords[0][3][1],
  bb_pts_coords[0][1][0],bb_pts_coords[0][1][1]
  );
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
function prepSrL8(image) {
  // Develop masks for unwanted pixels (fill, cloud, cloud shadow).
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);

  // Apply the scaling factors to the appropriate bands.
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
  // Replace original bands with scaled bands and apply masks.
  return image.addBands(scaled, null, true)
    .updateMask(qaMask).updateMask(saturationMask);
}
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
var l8_filtered = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filterDate(startDateTime, endDateTime)
  .map(prepSrL8)
  .select('SR.*','ST.*')
  .median();

var l8_clip = l8_filtered.clip(roi);

var l8_raw_filtered = ee.ImageCollection('LANDSAT/LC08/C02/T1')
  .filterDate(startDateTime, endDateTime);

// The asFloat parameter gives floating-point TOA output instead of
// the UINT8 outputs of the default simpleComposite().
var composite = ee.Algorithms.Landsat.simpleComposite({
  collection: l8_raw_filtered,
  asFloat: true
});

var l8_composite_clip = composite.clip(roi);
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// Display the cloud-free median composite.
var visParams = {
  //'SR_B6', 'SR_B5', 'SR_B3', '
  bands: ['ST_B10'],
  min: 0,
  max: 0.4
};

var thermalBands = l8_clip.select('ST_B10');
var opticalBands = l8_clip.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']);
var l8_composite_oli = l8_composite_clip.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7']);


/*var b1_dType = l8_composite_opticalBands.getInfo().bands[0].data_type;
var lst_dType = thermalBands.getInfo().bands[0].data_type;
print(b1_dType, lst_dType)*/

//var thermalBands = l8_clip.select('ST_B10').multiply(0.00341802).add(149.0);
var thermalBands_converted = thermalBands.expression(
  '(dn - 273.15) * 1.8 + 32.0', {
    'dn': thermalBands
  }
);

var thermalBands_converted = thermalBands_converted.rename('thermalBands_converted');

var LST_max = ee.Number(thermalBands_converted.reduceRegion({
  reducer: ee.Reducer.max(),
  geometry: roi,
  scale: 30,
  bestEffort: true,
  maxPixels: 1e7
}).values().get(0));


var LST_min = ee.Number(thermalBands_converted.reduceRegion({
  reducer: ee.Reducer.min(),
  geometry: roi,
  scale: 30,
  bestEffort: true,
  maxPixels: 1e7
}).values().get(0));

print('LST_Min', LST_min, LST_max); 

//Map.setCenter(mapCenter[0],mapCenter[1], 10);
//Map.addLayer(thermalBands.clip(roi), visParams, 'Cloud-free mosaic');
Map.addLayer(thermalBands_converted, {palette:['purple', 'blue', 'yellow', 'orange', 'red'], min: 75, max: 155}, 'thermalBands_converted');
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
var statistics = thermalBands_converted.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: roi,
  scale: 30,  // Adjust the scale as needed
  maxPixels: 1e9 // Specify maxPixels parameter as needed
});

var averageTemperature = statistics.get('thermalBands_converted');

// Print the average temperature
print('Average Temperature:', averageTemperature);

// Calculate the mode temperature
var mode = thermalBands_converted.reduce(ee.Reducer.mode());

// Get the mode temperature value
var modeTemperature = mode.get('thermalBands_converted');

// Print the mode temperature
print('Mode Temperature:', modeTemperature);
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//Export the l8_clip image
Export.image.toDrive({
  image: thermalBands_converted,
  description: place_key+'_L8_ARD_TIRS_'+year,
  folder: 'L8_ARD',
  region: roi,
  scale: 30,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: l8_composite_oli,
  description: place_key+'_L8_COMP_OLI_'+year,
  folder: 'L8_ARD',
  region: roi,
  scale: 30,
  fileFormat: 'GeoTIFF'
});
