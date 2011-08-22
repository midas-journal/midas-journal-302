#ifndef __itkQEMeshDiscreteCurvatureTensorEstimator_h
#define __itkQEMeshDiscreteCurvatureTensorEstimator_h

namespace itk
{
/**
 * \class QEMeshDiscreteCurvatureTensorEstimator
 * \brief
*/
template< class TInputMesh, class TOutputMesh >
class QEMeshDiscreteCurvatureTensorEstimator :
  public QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
{
public:
  typedef QEMeshDiscreteCurvatureTensorEstimator Self;
  typedef SmartPointer< Self > Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter Superclass;
  /** Run-time type information (and related methods).   */
  itkTypeMacro( QEMeshDiscreteCurvatureTensorEstimator,
    QuadEdgeMeshToQuadEdgeMeshFilter );
  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

protected:
  QEMeshDiscreteCurvatureTensorEstimator() : Superclass() {}
  ~QEMeshDiscreteCurvatureTensorEstimator() {}

  ///TODO to be implemented  
  virtual void GenerateData()
  {
    
  }
private:
  QEMeshDiscreteCurvatureTensorEstimator( const Self& );
  void operator = ( const Self& );
};
}
#endif
