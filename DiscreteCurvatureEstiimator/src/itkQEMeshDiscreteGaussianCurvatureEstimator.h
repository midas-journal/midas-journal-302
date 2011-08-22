#ifndef __itkQEMeshDiscreteGaussianCurvatureEstimator_h
#define __itkQEMeshDiscreteGaussianCurvatureEstimator_h

#include "itkQEMeshDiscreteCurvatureEstimator.h"
#include <vnl/vnl_math.h>

#include "itkTriangle.h"

namespace itk
{
/**
 * \class QEMeshDiscreteGaussianCurvatureEstimator
 * \brief see the following paper
 * title: Discrete Differential-Geometry Operators for Triangulated 2-Manifolds
 * authors: Mark Meyer, Mathieu Desbrun, Peter Schroder, Alan H. Barr
 * conference: VisMath '02
 * location: Berlin (Germany)
 * \author: Arnaud Gelas, Alexandre Gouaillard
*/
template< class TInputMesh, class TOutputMesh >
class QEMeshDiscreteGaussianCurvatureEstimator :
  public QEMeshDiscreteCurvatureEstimator< TInputMesh, TOutputMesh >
{
public:
  typedef QEMeshDiscreteGaussianCurvatureEstimator      Self;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;
  typedef QEMeshDiscreteCurvatureEstimator< TInputMesh, TOutputMesh >
    Superclass;

  typedef typename Superclass::InputMeshType    InputMeshType;
  typedef typename Superclass::InputMeshPointer InputMeshPointer;

  typedef typename Superclass::OutputMeshType   OutputMeshType;
  typedef typename Superclass::OutputMeshPointer OutputMeshPointer;
  typedef typename Superclass::OutputPointsContainerPointer 
    OutputPointsContainerPointer;
  typedef typename Superclass::OutputPointsContainerIterator
    OutputPointsContainerIterator;
  typedef typename Superclass::OutputPointType OutputPointType;
  typedef typename Superclass::OutputVectorType OutputVectorType;
  typedef typename Superclass::OutputCoordType OutputCoordType;
  typedef typename Superclass::OutputPointIdentifier OutputPointIdentifier;
  typedef typename Superclass::OutputCellIdentifier OutputCellIdentifier;
  typedef typename Superclass::OutputQEType OutputQEType;
  typedef typename Superclass::OutputMeshTraits OutputMeshTraits;
  typedef typename Superclass::OutputCurvatureType OutputCurvatureType;
  
  typedef typename Superclass::TriangleType TriangleType;
  
  /** Run-time type information (and related methods).   */
  itkTypeMacro( QEMeshDiscreteGaussianCurvatureEstimator,
    QEMeshDiscreteCurvatureEstimator );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

protected:
  QEMeshDiscreteGaussianCurvatureEstimator() : Superclass() {}
  ~QEMeshDiscreteGaussianCurvatureEstimator() {}

  virtual OutputCurvatureType EstimateCurvature( const OutputPointType& iP )
  {
    OutputMeshPointer output = this->GetOutput();

    OutputCurvatureType sum_theta;
    OutputCurvatureType area;

    OutputQEType* qe = iP.GetEdge( );
    OutputQEType* qe_it;
    OutputQEType* qe_it2;

    OutputPointType q0, q1;

    if( qe != 0 )
      {
      qe_it = qe;
      sum_theta = 0.;
      area = 0.;
      do
        {
        // cell_id = qe_it->GetLeft();
        qe_it2 = qe_it->GetOnext();
        q0 = output->GetPoint( qe_it->GetDestination() );
        q1 = output->GetPoint( qe_it2->GetDestination() );

        // Compute Angle;
        sum_theta += static_cast< OutputCurvatureType >(
              TriangleType::ComputeAngle( q0, iP, q1 ) );
        area += ComputeMixedArea( qe_it, qe_it2 );
        qe_it = qe_it2;
        } while( qe_it != qe );
    }
  return ( 2. * vnl_math::pi - sum_theta ) / area;
  }

private:
  QEMeshDiscreteGaussianCurvatureEstimator( const Self& );
  void operator = ( const Self& );
};
}
#endif
