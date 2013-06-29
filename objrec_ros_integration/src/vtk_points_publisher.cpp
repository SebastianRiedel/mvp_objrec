#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <stdint-gcc.h>

#include "ros/ros.h"
#include "sensor_msgs/PointCloud2.h"
#include "sensor_msgs/point_cloud_conversion.h"

// Helper function for raising an exception if a required parameter is not found
template <class T>
static void require_param(const ros::NodeHandle &nh, const std::string &param_name, T &var)
{
  if(!nh.getParam(param_name, var)) {
    ROS_FATAL_STREAM("Required parameter not found! Namespace: "<<nh.getNamespace()<<" Parameter: "<<param_name);
    throw ros::InvalidParameterException("Parameter not found!");
  }
}

int main(int argc, char *argv[])
{
  ros::init(argc,argv,"vtk_points_publisher");
  ros::NodeHandle nh;

  std::string inputFileName;
  double scale;
  require_param(nh,"vtk_scene_filename",inputFileName);
  require_param(nh,"vtk_scene_scale_to_meter",scale);

  ros::Publisher pub = nh.advertise<sensor_msgs::PointCloud2>("camera/depth_registered/points", 1);

  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(inputFileName.c_str());
  reader->Update();
  vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
  input = reader->GetOutput();

  sensor_msgs::PointCloud2 cloud;
  cloud.header.frame_id = "/camera_rgb_optical_frame";
  cloud.height = 1;
  cloud.width = input->GetNumberOfPoints();
  cloud.fields.resize(4);
  cloud.fields[0].name = "x";
  cloud.fields[0].offset = 0;
  cloud.fields[0].datatype = sensor_msgs::PointField::FLOAT32;
  cloud.fields[0].count = 1;
  cloud.fields[1].name = "y";
  cloud.fields[1].offset = 4;
  cloud.fields[1].datatype = sensor_msgs::PointField::FLOAT32;
  cloud.fields[1].count = 1;
  cloud.fields[2].name = "z";
  cloud.fields[2].offset = 8;
  cloud.fields[2].datatype = sensor_msgs::PointField::FLOAT32;
  cloud.fields[2].count = 1;
  cloud.fields[3].name = "rgb";
  cloud.fields[3].offset = 12;
  cloud.fields[3].datatype = sensor_msgs::PointField::FLOAT32;
  cloud.fields[3].count = 1;
  cloud.is_dense = false;
  cloud.is_bigendian = false;
  cloud.point_step = 16;
  cloud.row_step = cloud.point_step * cloud.width;
  cloud.data.resize(input->GetNumberOfPoints() * cloud.point_step);

  // copy points
  vtkPoints * pts = input->GetPoints();
  uint8_t *cloudPtr = &cloud.data[0];

  double pt[3];
  float pt_f[4];
  //cout << "#Points: " << pts->GetNumberOfPoints() << endl;
  for(int i=0; i<pts->GetNumberOfPoints(); ++i)
  {
    pts->GetPoint(i, pt);
    pt_f[0] = (float)pt[0]*scale;
    pt_f[1] = (float)pt[1]*scale;
    pt_f[2] = (float)pt[2]*scale;
    pt_f[3] = (float)(0xFFFFFFFF); // TODO: substitute with color from input vtk file
    memcpy(cloudPtr, &pt_f[0], cloud.point_step);
    cloudPtr += cloud.point_step;
  }

  ros::Rate loop_rate(30);
  ROS_INFO("Scaling factor is to: %f", scale);
  while(ros::ok())
  {
    pub.publish(cloud);
    ros::spinOnce();
    loop_rate.sleep();
  }

  return EXIT_SUCCESS;
}

