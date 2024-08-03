// Urban Growth Change Detection using Dynamic World Probability Bands
var expport = false;
/*var geometry = ee.Geometry.Polygon([[
  [77.43062052556523, 13.103764122826366],
  [77.43062052556523, 12.821384160047845],
  [77.7588370782996, 12.821384160047845],
  [77.7588370782996, 13.103764122826366]
]]);*/

var geometry = Filter_Polygon_B;

//Map.centerObject(geometry);

// Define the before and after time periods.
var beforeYear = 2018;
var afterYear = 2023;

// Create start and end dates for the before and after periods.
var beforeStart = ee.Date.fromYMD(beforeYear, 1 , 1);
var beforeEnd = beforeStart.advance(1, 'year');

var afterStart = ee.Date.fromYMD(afterYear, 1 , 1);
var afterEnd = afterStart.advance(1, 'year');

// Load the Dynamic World collection
var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')

// Filter the collection and select the 'built' band.
var dwFiltered = dw
  .filter(ee.Filter.bounds(geometry))
  .select('built');

// Create mean composites
var beforeDw = dwFiltered.filter(
  ee.Filter.date(beforeStart, beforeEnd)).mean();
  
var afterDw = dwFiltered.filter(
  ee.Filter.date(afterStart, afterEnd)).mean();


// Add Sentinel-2 Composites to verify the results.
var s2 = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
     .filterBounds(geometry)
     .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 35));

// Create a median composite from sentinel-2 images.
var beforeS2 = s2.filterDate(beforeStart, beforeEnd).median();
var afterS2 = s2.filterDate(afterStart, afterEnd).median();
  
// Visualize images
var s2VisParams = {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000};
Map.centerObject(geometry, 10);
Map.addLayer(beforeS2.clip(geometry), s2VisParams, 'Before S2');
Map.addLayer(afterS2.clip(geometry), s2VisParams, 'After S2');

// Select all pixels that have experienced large change
// in 'built' probbility
var builtChangeThreshold = 0.05; 
var newUrban = afterDw.subtract(beforeDw).gt(builtChangeThreshold);

var changeVisParams = {min: 0, max: 1, palette: ['white', 'red']};
Map.addLayer(newUrban.clip(geometry), changeVisParams, 'New Urban');

// Mask all pixels with 0 value using selfMask()
var newUrbanMasked = newUrban.selfMask();

Map.addLayer(
  newUrbanMasked.clip(geometry), changeVisParams, 'New Urban (Masked)');

// To ensure the masked values are set to NoData, 
// we cast the image to float and clip to geomery
var newUrbanMaskedExport = newUrbanMasked.toFloat().clip(geometry);

/*if (export === True){
  Export.image.toDrive({
    image: newUrbanMaskedExport,
    description: 'New_Urban_Areas_' + beforeYear + '_' + afterYear,
    folder: 'earthengine',
    fileNamePrefix: 'new_urban_areas_' + beforeYear + '_' + afterYear,
    region: geometry,
    scale: 10,
    maxPixels: 1e10
  });
};*/
