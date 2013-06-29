for now just some notes which eventually form a documentation

* loaded models (*.vtk and *.stl files) have to be in [mm] not [m]
* topic foreground_points publishes clipped scene, clipping in x,y,z direction in camera frame is set via ros params or dynamic reconfigure
* n_clouds_per_recognition is the maximum number of clouds stacked/added together before performing the recognition
* downsample voxel size is in mm

