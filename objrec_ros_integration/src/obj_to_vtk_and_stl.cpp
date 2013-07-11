#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkOBJReader.h>
#include <vtkSTLWriter.h>
#include <vtkPLYReader.h>
#include <vtkDecimatePro.h>
#include <vtkQuadricDecimation.h>
#include <vtkPolyDataWriter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>

/*
 * From http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataExtractNormals
 */
bool GetPointNormals(vtkPolyData* polydata)
{
  std::cout << "In GetPointNormals: " << polydata->GetNumberOfPoints() << std::endl;
  std::cout << "Looking for point normals..." << std::endl;

  // Count points
  vtkIdType numPoints = polydata->GetNumberOfPoints();
  std::cout << "There are " << numPoints << " points." << std::endl;

  // Count triangles
  vtkIdType numPolys = polydata->GetNumberOfPolys();
  std::cout << "There are " << numPolys << " polys." << std::endl;

  ////////////////////////////////////////////////////////////////
  // Double normals in an array
  vtkDoubleArray* normalDataDouble =
    vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetArray("Normals"));

  if(normalDataDouble)
    {
    int nc = normalDataDouble->GetNumberOfTuples();
    std::cout << "There are " << nc
            << " components in normalDataDouble" << std::endl;
    return true;
    }

  ////////////////////////////////////////////////////////////////
  // Double normals in an array
  vtkFloatArray* normalDataFloat =
    vtkFloatArray::SafeDownCast(polydata->GetPointData()->GetArray("Normals"));

  if(normalDataFloat)
    {
    int nc = normalDataFloat->GetNumberOfTuples();
    std::cout << "There are " << nc
            << " components in normalDataFloat" << std::endl;
    return true;
    }

  ////////////////////////////////////////////////////////////////
  // Point normals
  vtkDoubleArray* normalsDouble =
    vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetNormals());

  if(normalsDouble)
    {
    std::cout << "There are " << normalsDouble->GetNumberOfComponents()
              << " components in normalsDouble" << std::endl;
    return true;
    }

  ////////////////////////////////////////////////////////////////
  // Point normals
  vtkFloatArray* normalsFloat =
    vtkFloatArray::SafeDownCast(polydata->GetPointData()->GetNormals());

  if(normalsFloat)
    {
    std::cout << "There are " << normalsFloat->GetNumberOfComponents()
              << " components in normalsFloat" << std::endl;
    return true;
    }

  /////////////////////////////////////////////////////////////////////
  // Generic type point normals
  vtkDataArray* normalsGeneric = polydata->GetPointData()->GetNormals(); //works
  if(normalsGeneric)
    {
    std::cout << "There are " << normalsGeneric->GetNumberOfTuples()
              << " normals in normalsGeneric" << std::endl;

    double testDouble[3];
    normalsGeneric->GetTuple(0, testDouble);

    std::cout << "Double: " << testDouble[0] << " "
              << testDouble[1] << " " << testDouble[2] << std::endl;

    // Can't do this:
    /*
    float testFloat[3];
    normalsGeneric->GetTuple(0, testFloat);

    std::cout << "Float: " << testFloat[0] << " "
              << testFloat[1] << " " << testFloat[2] << std::endl;
    */
    return true;
    }


  // If the function has not yet quit, there were none of these types of normals
  std::cout << "Normals not found!" << std::endl;
  return false;

}

vtkSmartPointer<vtkPolyData> decimate(vtkSmartPointer<vtkPolyData> input, double targetReduction)
{
  vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
  decimate->SetInput(input);
  decimate->SetTargetReduction(targetReduction);
  decimate->SetPreserveTopology(1);
  decimate->Update();

  vtkSmartPointer<vtkPolyData> decimated = decimate->GetOutput();
  return decimated;
}

bool loadModel(std::string filename, vtkSmartPointer<vtkPolyData> &model)
{
  vtkSmartPointer<vtkPolyData> inputMesh = vtkSmartPointer<vtkPolyData>::New();
  int sep = filename.find_last_of('.');
  if(filename.substr(sep) == ".vtk")
  {
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    inputMesh = reader->GetOutput();
  }
  else if(filename.substr(sep) == ".obj")
  {
    vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    inputMesh = reader->GetOutput();
  }
  else if(filename.substr(sep) == ".ply")
  {
    vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    inputMesh = reader->GetOutput();
  }
  else
  {
    std::cerr << "Can't load model of type: " << filename.substr(sep) << std::endl;
    return false;
  }
  model = inputMesh;
  return true;
}

int main(int argc, char *argv[])
{
  if(argc != 6)
  {
      std::cout << "Usage: obj_to_vtk_and_stl InputObjFile OutputStlFile OutputVtkFile maxNPointsMeshShouldHaveAfterDecimation ScaleToConvertObjMeshToMM" << std::endl;
      return EXIT_SUCCESS;
  }
  std::string inputFileName = argv[1];
  std::string outputFileNameStl = argv[2];
  std::string outputFileNameVtk = argv[3];
  int maxNPoints = atoi(argv[4]);
  float scale = atof(argv[5]);

  std::cout << "Reading " << inputFileName << std::endl;
  vtkSmartPointer<vtkPolyData> input;
  loadModel(inputFileName, input);

  std::cout << "Input has normals: " << GetPointNormals(input) << std::endl;
  if(!GetPointNormals(input))
  {
    std::cout << "WARNING: Input mesh, does not have normals! Resulting vtk file will not have normals." << std::endl;
  }

  int maxIterations = 25;
  int curIteration = 0;
  std::cout << "Trying to decimate mesh to maximal  " << maxNPoints << "points within maximal " << maxIterations << " iterations." << std::endl;
  std::cout << ".. " << curIteration << ": NPoints " << input->GetNumberOfPoints() << std::endl;
  while(input->GetNumberOfPoints() > maxNPoints && ++curIteration < maxIterations)
  {
    {
    vtkSmartPointer<vtkPolyData> temp = decimate(input, 0.5);
    input = temp;
    }
    std::cout << ".. " << curIteration << ": NPoints " << input->GetNumberOfPoints() << std::endl;
  }

  std::cout << "Applying scaling by factor or " << scale << std::endl;
  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->Scale(scale,scale,scale);
  vtkSmartPointer<vtkTransformPolyDataFilter> scaleFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  scaleFilter->SetInput(input);
  scaleFilter->SetTransform(transform);
  scaleFilter->Update();
  vtkSmartPointer<vtkPolyData> transformed = scaleFilter->GetOutput();


  std::cout << "Saving stl..." << std::endl;
  vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
  writer->SetFileName(outputFileNameStl.c_str());
  writer->SetInput(transformed);
  writer->SetFileTypeToBinary();
  writer->Update();

  std::cout << "Saving vtk..." << std::endl;
  vtkSmartPointer<vtkPolyDataWriter> writerVtk = vtkSmartPointer<vtkPolyDataWriter>::New();
  writerVtk->SetFileTypeToASCII();
  writerVtk->SetFileName(outputFileNameVtk.c_str());
  writerVtk->SetInput(transformed);
  writerVtk->Update();

  if(!GetPointNormals(transformed))
  {
      std::cout << "WARNING: Output meshes do not contain normals (likely because input file did not have normals), detection won't work!!!" << std::endl;
  }

  return EXIT_SUCCESS;
}

