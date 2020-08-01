# 3d-maps-Guwahati-Satellite-Imagery-Landsat-8

![Demo1](https://github.com/Hwoabam/3d-maps-Guwahati-Satellite-Imagery-Landsat-8/blob/master/Media/Animation/GIF1..gif)

Various packages are required for the 3D map generation which are installed and loaded into the program. These packages avail the functions such as Plotting of points, shading and adjusting the colors, conversion of GIS data, validating GDAL operations,etc. 
```{r}
library(rayshader)
library(magick)
library(sp)
library(raster)
library(scales)
library(rgdal)
library(rgeos)
```
The GIS data for elevation is downloaded from Derek Watkin's -"30 meter SRTM tile downloader". The downloaded SRTM hgt data is converted to a matrix. The result of the plot is an elevation intensity profile of the tile.
```{r fig1, fig.height = 12, fig.width = 8, align= "center"}
windowsize=c(12,8)

guwahati_elevation = raster("D:/Assam maps/New Folder/N26E091.hgt")
height_shade(raster_to_matrix(guwahati_elevation)) %>%
  plot_map()
```
![Elevation Heat Map](https://github.com/Hwoabam/3d-maps-Guwahati-Satellite-Imagery-Landsat-8/blob/master/Media/Plots/Elevation_heatmap.png)

Now the satellite data is required to be downloaded from the USGS Earth explorer or any other suitable sources(Satellite 8 Data is specifically required for this purpose). The B2, B3 and B4 implies to blue, green and red color bands. So these TIF files are loaded and then stacked for band combinations and a satellite image is plotted 
```{r fig2, fig.height = 12, fig.width = 8, align= "center"}
guwahati_r = raster("D:/Assam maps/LC08_L1TP_137042_20200330_20200409_01_T1_B4.TIF")
guwahati_g = raster("D:/Assam maps/LC08_L1TP_137042_20200330_20200409_01_T1_B3.TIF")
guwahati_b = raster("D:/Assam maps/LC08_L1TP_137042_20200330_20200409_01_T1_B2.TIF")

guwahati_rgb = stack(guwahati_r, guwahati_g, guwahati_b)
plotRGB(guwahati_rgb, scale=255^2)
```
![Truecolor plot](https://github.com/Hwoabam/3d-maps-Guwahati-Satellite-Imagery-Landsat-8/blob/master/Media/Plots/RGB.png)

The above plot lacks brightness which is increased by implementing Gamma correction
```{r fig3, fig.height = 12, fig.width = 8, align= "center"}
guwahati_rgb_corrected = sqrt(stack(guwahati_r, guwahati_g, guwahati_b))
plotRGB(guwahati_rgb_corrected)
```
![Gamma correction plot](https://github.com/Hwoabam/3d-maps-Guwahati-Satellite-Imagery-Landsat-8/blob/master/Media/Plots/landsat_8_gamma_corrected.png)

The coordinate reference of the elevation data matched to the imagery data
```{r}
crs(guwahati_r)
crs(guwahati_elevation)
```
Transforming data from Long/lat to UTM using 'bilinear' method of interpolation. since here elevation is a continuous variable 
```{r}
guwahati_elevation_utm = projectRaster(guwahati_elevation, crs = crs(guwahati_r), method = "bilinear")
crs(guwahati_elevation_utm)
```
The desired area whose map is to be plot is cropped in a rectangular shape by defining the co-ordinates of the bottom left and top right, diagonally opposite points. The longitude and latitude of the two points are found using "Google Maps" and then defined.
```{r}
bot_lefty=91.529469
bot_leftx=26.065919
top_righty=91.870045
top_rightx=26.305610
distl2l= top_righty-bot_lefty  
scl=1-(20/(distl2l*111))
bottom_left = c(bot_lefty, bot_leftx)
top_right   = c(top_righty, top_rightx)
extent_latlong = SpatialPoints(rbind(bottom_left, top_right), proj4string=sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
extent_utm = sp::spTransform(extent_latlong, raster::crs(guwahati_elevation_utm))
e = raster::extent(extent_utm)
e
scl
```
Now the Cropped area is to be transformed into a 3 layered RGB array and the elevations need to be transposed on the array. The aperm() function is used for this purpose.Then, the elevation data is converted into Base R matrix
```{r fig4, fig.height = 12, fig.width = 8, align= "center"}
guwahati_rgb_cropped = crop(guwahati_rgb_corrected, e)
elevation_cropped = crop(guwahati_elevation_utm, e)
names(guwahati_rgb_cropped) = c("r","g","b")
guwahati_r_cropped = raster_to_matrix(guwahati_rgb_cropped$r)
guwahati_g_cropped = raster_to_matrix(guwahati_rgb_cropped$g)
guwahati_b_cropped = raster_to_matrix(guwahati_rgb_cropped$b)
guwahatiel_matrix = raster_to_matrix(elevation_cropped)
guwahati_rgb_array = array(0,dim=c(nrow(guwahati_r_cropped),ncol(guwahati_r_cropped),3))
guwahati_rgb_array[,,1] = guwahati_r_cropped/255 #Red layer
guwahati_rgb_array[,,2] = guwahati_g_cropped/255 #Blue layer
guwahati_rgb_array[,,3] = guwahati_b_cropped/255 #Green layer
guwahati_rgb_array = aperm(guwahati_rgb_array, c(2,1,3))
plot_map(guwahati_rgb_array)
```
![Cropped matrix plot](https://github.com/Hwoabam/3d-maps-Guwahati-Satellite-Imagery-Landsat-8/blob/master/Media/Plots/landsat_8_crop.png)

The contrast is then adjusted
```{r fig5, fig.height = 12, fig.width = 8, align= "center"}
guwahati_rgb_contrast = rescale(guwahati_rgb_array,to=c(0,1))
plot_map(guwahati_rgb_contrast)
```
![Final contrast corrected plot](https://github.com/Hwoabam/3d-maps-Guwahati-Satellite-Imagery-Landsat-8/blob/master/Media/Plots/landsat_8_crop_contrast_corrected.png)

Plotting of 3d model in an rgl window. Defining viewpoint parameters in addition to Zscale which is set at 7.5 for more distinct topography. Background color, Shadow color are defined too with default color codes. Compass rendering is done with its position to the West of the map.
```{r fig6, fig.height = 15, fig.width = 10, align= "center"}
plot_3d(guwahati_rgb_contrast, guwahatiel_matrix, windowsize = c(1100,900), zscale = 7.5, shadowdepth = -150,
        zoom=0.5, phi=45,theta=-45,fov=70, background = "#F2E1D0", shadowcolor = "#523E2B")
render_compass(position = "W", )
render_snapshot(title_text = "Guwahati City | Imagery: Landsat 8 | DEM: 30m SRTM", title_bar_color = "#1f5214", title_color = "white", title_bar_alpha = 1)
```
![3D plot with DEM and Satellite Imagery ](https://github.com/Hwoabam/3d-maps-Guwahati-Satellite-Imagery-Landsat-8/blob/master/Media/Snapshots/snap1.png)

Conversion of the map to a video footage of its planar rotation using snapshots taken at 1440 intervals during the rotation at 60 fps with the help of ffmpeg
```{r}
angles= seq(0,360,length.out = 1441)[-1]
for(i in 1:1440) {
  render_camera(theta=-45+angles[i])
  render_snapshot(filename =
  sprintf("guwahati_sat%i.png", i), 
                  title_text = "Guwahati_Satellite_Imagery| DEM: 30m SRTM",
                  title_bar_color = "#1f5214", title_color = "white", title_bar_alpha = 1)
}
rgl::rgl.close()
system("ffmpeg -framerate 60 -i guwahati_sat%d.png -vcodec libx264 -an Guwahati_satellite.mp4 ")
```
![Demo1](https://github.com/Hwoabam/3d-maps-Guwahati-Satellite-Imagery-Landsat-8/blob/master/Media/Animation/GIF1..gif)

The geographical plot with 3d visualization in an rgl window. Defining viewpoint parameters in addition to Zscale which is set at 7.5 for more distinct topography. Background color, Shadow color are defined too with default color codes. Compass rendering is done with its position to the West of the map. Scalebar rendering is done with markings at 0,5, 10 (km) upto the caluclated length from the position 1.
```{r fig7, fig.height = 15, fig.width = 10, align= "center"}
guwahatiel_matrix %>%
  sphere_shade(texture = "imhof1") %>%
  add_water(detect_water(guwahatiel_matrix), color = "imhof4") %>%
  plot_3d( guwahatiel_matrix, windowsize = c(1100,900), zscale = 7.5, shadowdepth = -150,
           zoom=0.5, phi=45,theta=-45,fov=70, background = "#F2E1D0", shadowcolor = "#523E2B")
render_scalebar(limits=c(20,10,0),label_unit = "km",position = "S", y=50,scale_length = c(scl,1))
render_compass(position = "W" )
render_snapshot(title_text = "Guwahati_Geographical|DEM: 30m SRTM",title_bar_color = "#1f5214", title_color = "white", title_bar_alpha = 1)
```
![3D plot with DEM and Satellite Imagery ](https://github.com/Hwoabam/3d-maps-Guwahati-Satellite-Imagery-Landsat-8/blob/master/Media/Snapshots/snap2.png)

![Demo2](https://github.com/Hwoabam/3d-maps-Guwahati-Satellite-Imagery-Landsat-8/blob/master/Media/Animation/GIF2.gif)

Conversion of the map to a video footage of its planar rotation using snapshots taken at 1440 intervals during the rotation at 60 fps with the help of ffmpeg
```{r}
angles= seq(0,360,length.out = 1441)[-1]
for(i in 1:1440) {
render_camera(theta=-45+angles[i])
  render_snapshot(filename = sprintf("guwahati_geo%i.png", i), 
                  title_text = "Guwahati_geographical| DEM: 30m SRTM",
                  title_bar_color = "#1f5214", title_color = "white", title_bar_alpha = 1)
}
rgl::rgl.close()
system("ffmpeg -framerate 60 -i guwahati_geo%d.png -vcodec libx264 -an Guwahati_geographical_1.mp4 ")
```
