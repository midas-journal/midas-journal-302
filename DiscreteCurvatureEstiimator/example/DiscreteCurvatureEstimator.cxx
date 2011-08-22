#include <itkQuadEdgeMesh.h>
#include <itkVTKPolyDataReader.h>

#include "itkQuadEdgeMeshExtendedTraits.h"
#include "itkQEMeshDiscreteGaussianCurvatureEstimator.h"
#include "itkQEMeshDiscreteMeanCurvatureEstimator.h"
#include "itkQEMeshDiscreteMinCurvatureEstimator.h"
#include "itkQEMeshDiscreteMaxCurvatureEstimator.h"
#include "itkQEMeshScalarDataVTKPolyDataWriter.h"

using namespace itk;

int main( int argc, char** argv )
{
  if( argc < 3 )
    {
    std::cout <<"*** Discrete Curvature Estimator ***" <<std::endl;
    std::cout <<"This example requires at least one argument:" <<std::endl;
    std::cout <<" 1- FileName" <<std::endl;
    std::cout <<" 2- CurvatureType" <<std::endl;
    std::cout <<"   * 0: Gaussian" <<std::endl;
    std::cout <<"   * 1: mean" <<std::endl;
    std::cout <<"   * 2: min" <<std::endl;
    std::cout <<"   * 3: max" <<std::endl;
    std::cout <<"************************************" <<std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef double CoordType;

  typedef QuadEdgeMeshExtendedTraits <
    CoordType,
    Dimension,
    2,
    CoordType,
    CoordType,
    CoordType,
    bool,
    bool > Traits;

  typedef QuadEdgeMesh< CoordType, Dimension, Traits > MeshType;
  typedef QEMeshDiscreteGaussianCurvatureEstimator<MeshType,MeshType>
    GaussianCurvatureFilterType;
  typedef QEMeshDiscreteMeanCurvatureEstimator<MeshType,MeshType>
    MeanCurvatureFilterType;
  typedef QEMeshDiscreteMinCurvatureEstimator<MeshType,MeshType>
    MinCurvatureFilterType;
  typedef QEMeshDiscreteMaxCurvatureEstimator<MeshType,MeshType>
    MaxCurvatureFilterType;

  typedef VTKPolyDataReader< MeshType > ReaderType;

  ReaderType::Pointer reader = ReaderType::New( );
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update( );
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  MeshType::Pointer mesh = reader->GetOutput();

  GaussianCurvatureFilterType::Pointer
    gaussian_curvature = GaussianCurvatureFilterType::New();

  MeanCurvatureFilterType::Pointer
    mean_curvature = MeanCurvatureFilterType::New();

  MinCurvatureFilterType::Pointer
    min_curvature = MinCurvatureFilterType::New();

  MaxCurvatureFilterType::Pointer max_curvature =
        MaxCurvatureFilterType::New();

  MeshType::Pointer output;
  std::string output_filename;

  switch( atoi( argv[2] ) )
  {
    case 0:
      gaussian_curvature->SetInput( mesh );
      gaussian_curvature->Update();

      output = gaussian_curvature->GetOutput();
      output_filename = "gaussian_curvature.vtk";
      break;
    case 1:
      mean_curvature->SetInput( mesh );
      mean_curvature->Update();

      output = mean_curvature->GetOutput();
      output_filename = "mean_curvature.vtk";
      break;
    case 2:
      min_curvature->SetInput( mesh );
      min_curvature->Update();

      output = min_curvature->GetOutput();
      output_filename = "min_curvature.vtk";
      break;
    case 3:
      max_curvature->SetInput( mesh );
      max_curvature->Update();

      output = max_curvature->GetOutput();
      output_filename = "max_curvature.vtk";
      break;
    default:
      std::cout <<"The second parameter should be in between 0 and 3"
        <<std::endl;
      return EXIT_FAILURE;
  }

  typedef QEMeshScalarDataVTKPolyDataWriter< MeshType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( output );
  writer->SetFileName( output_filename );
  writer->Update();

  return EXIT_SUCCESS;
}
