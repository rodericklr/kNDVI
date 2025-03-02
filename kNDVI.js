//kNDVI was calculated using Javascript based on GEE platform
//kNDVI uses the idea of kernel function,seeWang X, Biederman J A, Knowles J F, et al. Satellite solar-induced chlorophyll fluorescence and near-infrared reflectance capture complementary aspects of dryland vegetation productivity dynamics[J].Remote sensing of environment, 2022, 270: 112858.
//Because the research area is too large, the research area is divided into 4 small areas by segmentation function,and kNDVI calculation of each subarea is carried out successively.
//In order to make full use of the available effective observation pixels, kNDVI in 2000 was calculated based on LandSAT 5/7, and calculate kNDVI in 2023 based on Landsat8/9.
//Although the research area was divided, the results could not be displayed in GEE due to memory limitations. Therefore, this demo randomly selects a small area in the research area for code running and display

////////////////////////////////////////////////Step 1:The large study area is divided into several subregions/////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define the coordinates of the study area
var coordinates = [
  [24.818431891686092,-6.909813631787851],
  [27.976963112145977,-6.909813631787851],
  [27.976963112145977,-4.064109347062295],
  [24.818431891686092,-4.064109347062295],
  [24.818431891686092,-6.909813631787851] // 闭合矩形
];

//Create a rectangular ee.Geometry
var roi = ee.Geometry.Polygon([coordinates]);

// Set the style
var style = {color: 'black', fillColor: '00000000', width: 2};

// Add the rectangle to the map
Map.centerObject(roi, 5);
//Map.addLayer(roi, style, 'TP');

// Define a function to slice the study area and divide the study into n*m blocks
function gridSplit (roi, n, m) {
var bounds = roi.bounds().getInfo().coordinates[0];
//In order to ensure that the edge pixels are included in the range 
//and improve the fault tolerance, the range is expanded by 0.1 to all sides
var xmin = bounds[0][0]-0.1;
var xmax = bounds[1][0]+0.1;
var ymin = bounds[0][1]-0.1;
var ymax = bounds[2][1]+0.1;
//print(xmin,ymin,xmax,ymax);
// Calculates the unit length of the block
var dx = (xmax-xmin)/n;
var dy = (ymax-ymin)/m;
//Loop generates blocks
var xx = ee.List.sequence(xmin, xmax, dx);
var yy = ee.List.sequence(ymin, ymax, dy);
//To simplify the code, some blocks are generated that are beyond the scope of the study area
var rects = xx.map(function(i){
  return yy.map(function(j){
  var x1 = ee.Number(i);
  var x2 = ee.Number(i).add(ee.Number(dx));
  var y1 = ee.Number(j);
  var y2 = ee.Number(j).add(ee.Number(dy));
  var coords = ee.List([x1, y1, x2, y2]);
  var rect = ee.Algorithms.GeometryConstructors.Rectangle(coords);
  return ee.Feature(rect);
  });
}).flatten();
//Screen out the blocks intersecting the study area
var rects_col = ee.FeatureCollection(rects).filterBounds(roi);
//Calculate how many blocks there are in total
var GridNum =rects_col.size().getInfo();
//print('GridNum: ', GridNum);
//Assign a number to each block
var idList=ee.List.sequence(0, GridNum-1);
var grid = ee.FeatureCollection(idList.map(function(i) {
  return ee.Feature(rects_col.toList(rects_col.size()).get(i)).set("grid_id",ee.Number(i).add(1));
}));
return grid;
}
//The study was divided into 2*2 blocks, which were divided into 4 small areas
var subroi = gridSplit(roi, 2, 2);
//Visualize all blocks
Map.addLayer(subroi,{}, 'subroi_all');
//print(subroi,'subroi_list');

var roi2=subroi.filterMetadata('grid_id','equals', 1)//Modify this number to select a boundary range for different partitions(1,2,3,4)
Map.addLayer(roi2,{}, 'subroi_1')

/////////////////////////////////////////////Step2:2000 year kNDVI caclulation//////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// L57 cloud mask
function maskL57sr(image) {
      // Bit 0 - Fill
      // Bit 1 - Dilated Cloud
      // Bit 2 - Unused
      // Bit 3 - Cloud
      // Bit 4 - Cloud Shadow
      var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
      var saturationMask = image.select('QA_RADSAT').eq(0);
    
      // Apply the scaling factors to the appropriate bands.
      var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
      var thermalBand = image.select('ST_B6').multiply(0.00341802).add(149.0);
    
      // Replace the original bands with the scaled ones and apply the masks.
      return image.addBands(opticalBands, null, true)
          .addBands(thermalBand, null, true)
          .updateMask(qaMask)
          .updateMask(saturationMask);
    }
    

//Landsat57 kNDVI
var apply_sigma2_l57 = function(img){
    // Check if 'SR_B5' band is present and not null
  
  var red = img.select('SR_B3')
  var nir = img.select('SR_B4')
  // D^2
    //var D2 = nir.subtract(red).pow(2);
    var sigma2 = nir.add(red).multiply(0.5).reduceRegion({
              reducer:ee.Reducer.mean(), 
              geometry:roi2,
                scale: 500,
                crs: 'EPSG:4326',
                maxPixels:10e16,
                tileScale:16}).values().get(0);
  //var D2 = nir.subtract(red).pow(2);

    return img.set('sigma2',sigma2)
}

var apply_kNDVI_l57 = function(img){
    // Check if 'SR_B5' band is present and not null
  
  var red = img.select('SR_B3')
  var nir = img.select('SR_B4')
  // D^2
    var sigma2 = ee.Number(img.get('sigma2'))
    var D2 = nir.subtract(red).pow(2);
  //print('sigma', sigma2);
  sigma2 = ee.Number(sigma2).pow(2.0).multiply(2.0);
  // print('Estimated by 0.5*(nir+red)', sigma2);
  
  // k := exp(-D^2/sigma2)
  var k = D2.divide(sigma2).multiply(-1.0).exp();
  
  // kNDVI = (1-k)/(1+k)
  var kndvi = (ee.Image.constant(1).subtract(k))
      .divide(ee.Image.constant(1).add(k));
  return img.addBands(kndvi.rename('kNDVI'))
}


var l5Col = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
              .filterBounds(roi2)
              .filterDate("2000-04-01","2000-10-31")
              .map(maskL57sr)
              .map(apply_sigma2_l57)
              
var l5Col_filter = l5Col.filter(ee.Filter.neq('sigma2', null))    
                        .map(apply_kNDVI_l57)
                        .select('kNDVI')
                        .filterDate("2000-04-01","2000-10-31")
                        .max()
                       .clip(roi2) 

var l7Col = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
              .filterBounds(roi2)
              .filterDate("2000-04-01","2000-10-31")//Change the year to 2023 when running kNDVI 2023
              .map(maskL57sr)
              .map(apply_sigma2_l57)
              
var l7Col_filter = l7Col.filter(ee.Filter.neq('sigma2', null))    
                        .map(apply_kNDVI_l57)
                        .select('kNDVI')
                        .filterDate("2000-04-01","2000-10-31")
                        .max()
                        .clip(roi2)   
                      
var Col_NDVI2000 = ee.ImageCollection([l5Col_filter,l7Col_filter]).max();//for 2021

//Show kNDVI     
var ndviVis = {
min: -1,  // 颜色映射的最小值
max: 1,  // 颜色映射的最大值
palette: [
'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
'66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
'012E01', '011D01', '011301'
],
};
Map.addLayer(Col_NDVI2000,ndviVis,"Col_NDVI2000")

//Export kNDVI to drive
Export.image.toAsset({
  image: Col_NDVI2000.multiply(1000000).toUint32(),//Store in toUint32 format to reduce memory 
  description: 'kNDVI2000',
  assetId: 'kNDVI',
  scale: 30,
  maxPixels:1e13,
  region: roi2,
  });

/////////////////////////////////Step3: 2023 year kNDVI //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//L89 cloud mask
function applyScaleFactorsL89(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}
function cloudmaskL89(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 4);
  var cloudsBitMask = (1 << 3);
  // Get the pixel QA band.
  var qa = image.select('QA_PIXEL');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}


//Landsat89 kNDVI
var apply_sigma2_l89 = function(img){
    // Check if 'SR_B5' band is present and not null
  
  var red = img.select('SR_B4')
  var nir = img.select('SR_B5')
  // D^2
    var D2 = nir.subtract(red).pow(2);
    var sigma2 = nir.add(red).multiply(0.5).reduceRegion({
              reducer:ee.Reducer.mean(), 
              geometry:roi2,
                scale: 500,
                crs: 'EPSG:4326',
                maxPixels:10e16,
                tileScale:16}).values().get(0);
  var D2 = nir.subtract(red).pow(2);

    return img.set('sigma2',sigma2)
}

var apply_kNDVI_l89 = function(img){
    // Check if 'SR_B5' band is present and not null
  
  var red = img.select('SR_B4')
  var nir = img.select('SR_B5')
  // D^2
    var sigma2 = ee.Number(img.get('sigma2'))
    var D2 = nir.subtract(red).pow(2);
  //print('sigma', sigma2);
  sigma2 = ee.Number(sigma2).pow(2.0).multiply(2.0);
  // print('Estimated by 0.5*(nir+red)', sigma2);
  
  // k := exp(-D^2/sigma2)
  var k = D2.divide(sigma2).multiply(-1.0).exp();
  
  // kNDVI = (1-k)/(1+k)
  var kndvi = (ee.Image.constant(1).subtract(k))
      .divide(ee.Image.constant(1).add(k));
  return img.addBands(kndvi.rename('kNDVI'))
}


//Landsat 8 kNDVI caculation
var l8Col = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
              .filterBounds(roi2)
              .filterDate("2023-04-01", "2023-10-31")
              .map(applyScaleFactorsL89)
              .map(cloudmaskL89)
              .map(apply_sigma2_l89)
              
var l8Col_filter = l8Col.filter(ee.Filter.neq('sigma2', null))    
                                  .map(apply_kNDVI_l89)
                                  .select('kNDVI')
                                  .filterDate("2023-04-01","2023-10-31")
                                  .max()
                                  .clip(roi2)
 
//Landsat9 kNDVI caculation
var l9Col = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
                  .filterBounds(roi2)
                  .filterDate("2023-04-01", "2023-10-31")
                  .map(applyScaleFactorsL89)
                  .map(cloudmaskL89)
                  .map(apply_sigma2_l89);
                  
var l9Col_filter = l9Col.filter(ee.Filter.neq('sigma2', null))    
                                  .map(apply_kNDVI_l89)
                                  .select('kNDVI')
                                  .filterDate("2023-04-01", "2023-10-31")
                                  .max()
                                  .clip(roi2)

var Col_NDVI2023 = ee.ImageCollection([l8Col_filter,l9Col_filter]).max();//for 2023 
//Export kNDVI to drive
Export.image.toAsset({
  image: Col_NDVI2023.multiply(1000000).toUint32(),
  description: 'kNDVI2023',
  assetId: 'kNDVI',
  scale: 30,
  maxPixels:1e13,
  region: roi2,
  });
Map.addLayer(Col_NDVI2023,ndviVis,"Col_NDVI2023")  
