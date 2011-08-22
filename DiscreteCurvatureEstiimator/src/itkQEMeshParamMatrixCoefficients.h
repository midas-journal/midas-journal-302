#ifndef __itkQEMeshParamMatrixCoefficients_h
#define __itkQEMeshParamMatrixCoefficients_h

#include "itkQuadEdgeMesh.h"
#include "itkCross.h"
#include "itkTriangle.h"
#include <vnl/vnl_math.h>

namespace itk
{
/**
* \brief Superclass for all the matrix coefficients computation classes.
* \note  Belongs to the parameterisation package.
*/
template< typename TInputMesh >
class MatrixCoefficients
{
public:
    typedef TInputMesh InputMeshType;
    typedef typename InputMeshType::CoordRepType    InputCoordRepType;
    typedef typename InputMeshType::QEType          InputQEType;

    MatrixCoefficients( ){}
    virtual ~MatrixCoefficients( ) {}

    virtual InputCoordRepType operator ( )
        ( InputMeshType* iMesh, InputQEType* iEdge ) const = 0;
};

/**
* \brief Compute a matrix filled by 1s wherever two vertices are connected
*        by an edge.
* \note  Belongs to the parameterisation package.
* \note  See paper:
*/
template< typename TInputMesh >
class OnesMatrixCoefficients : public MatrixCoefficients< TInputMesh >
{
public:
    typedef MatrixCoefficients< TInputMesh > Superclass;

    typedef TInputMesh InputMeshType;
    typedef typename InputMeshType::CoordRepType    InputCoordRepType;
    typedef typename InputMeshType::QEType          InputQEType;

    OnesMatrixCoefficients( ) : Superclass( )
    { };

    /**
    * \param[in] iMesh
    * \param[in] iEdge
    * \return \f$ 1 \f$
    */
    InputCoordRepType operator ( ) ( InputMeshType* iMesh,
                                    InputQEType* iEdge ) const
    {
        return 1.;
    };
};

/**
* \brief Compute a matrix filed with the inverse of the euclidian distance
*        wherever two vertices are connected by an edge.
* \note  Belongs to the parameterisation package.
* \note  See paper: ...
*/
template< typename TInputMesh >
class InverseEuclideanDistanceMatrixCoefficients :
    public MatrixCoefficients< TInputMesh >
{
public:
    typedef MatrixCoefficients< TInputMesh > Superclass;

    typedef TInputMesh InputMeshType;
    typedef typename InputMeshType::CoordRepType    InputCoordRepType;
    typedef typename InputMeshType::PointType       InputPointType;
    typedef typename InputMeshType::PointIdentifier InputPointIdentifier;
    typedef typename InputMeshType::QEType          InputQEType;
    typedef typename InputMeshType::VectorType      InputVectorType;

    InverseEuclideanDistanceMatrixCoefficients( ) :
        Superclass( )
    { }

    /**
    * \param[in] iMesh
    * \param[in] iEdge
    * \return \f$ \frac{1}{\|\boldsymbol{p1} - \boldsymbol{p2} \|} \f$
     */
    InputCoordRepType operator ( ) ( InputMeshType* iMesh,
                                InputQEType* iEdge ) const
    {
        InputPointIdentifier id1 = iEdge->GetOrigin( );
        InputPointIdentifier id2 = iEdge->GetDestination( );

        InputPointType pt1 = iMesh->GetPoint( id1 );
        InputPointType pt2 = iMesh->GetPoint( id2 );

        InputCoordRepType oValue = 1. / pt1.EuclideanDistanceTo( pt2 );

        return oValue;
    }
};

//    /**
//    * \brief Compute the cross product of two vectors of dimension 3,
//    *        independently of the type of the values of vector's elements.
//    * \note  Belongs to the parameterisation package.
//    */
// template< typename TVector >
// class Cross
// {
// public:
//     typedef TVector VectorType;
//     typedef typename VectorType::ValueType ValueType;
// 
//     itkStaticConstMacro( Dimension, unsigned int, VectorType::Dimension );
// 
//     /**
//      * \param[in] iU
//      * \param[in] iV
//      * \return \f$ \boldsymbol{iU} \cdot \boldsymbol{iV} \f$
//      */
//     VectorType operator ( ) ( const VectorType& iU,
//                                 const VectorType& iV ) const
//     {
//         assert( Dimension == 3 );
// 
//         VectorType oCross;
// 
//         // *** Arnaud : InputVDimension == 3
//         oCross[0] = iU[1] * iV[2] - iV[1] * iU[2];
//         oCross[1] = iV[0] * iU[2] - iU[0] * iV[2];
//         oCross[2] = iU[0] * iV[1] - iV[0] * iU[1];
// 
//         return oCross;
//     }
// };

//    /**
//    * \brief Compute the Cotangent of an angle given by three points,
//    *        independently of the type of coordinates (int, float, double,...).
//    * \note  Used by the parameterisation package.
//    */
// template< typename TPoint >
// class Cotangent
// {
// public:
//     typedef TPoint PointType;
//     typedef typename PointType::CoordRepType  CoordRepType;
// 
//     itkStaticConstMacro( PointDimension, unsigned int,
//         PointType::PointDimension );
// 
//     typedef typename PointType::VectorType VectorType;
// 
//     /**
//     * \param[in] iPt1
//     * \param[in] iPt2
//     * \param[in] iPt3
//     * \return \f$ \text{cot}\left(\hat{iPt1\ iPt2\ iPt3}\right) \f$
//     */
//     CoordRepType operator ( ) ( const PointType& iPt1,
//                                 const PointType& iPt2,
//                                 const PointType& iPt3 ) const
//     {
//         VectorType v21 = iPt1 - iPt2;
//         CoordRepType squared_norm21 = v21.GetSquaredNorm( );
// 
//         VectorType v23 = iPt3 - iPt2;
//         CoordRepType squared_norm23 = v23.GetSquaredNorm( );
// 
//         CoordRepType squared_cos_theta = v21 * v23;
//         squared_cos_theta *= squared_cos_theta / ( squared_norm21 *
// squared_norm23 );
//         CoordRepType sin_theta = vcl_sqrt( 1. - squared_cos_theta );
// 
//         CoordRepType oValue = vcl_sqrt( squared_cos_theta ) / sin_theta;
// 
//         return oValue;
//     };
// };

/**
    * \brief Compute a matrix filed by Conformal Coefficients of the edge
    *        wherever two vertices are connected by an edge.
    * \note  Belongs to the parameterisation package.
    * \note  See paper ...
    */
template< typename TInputMesh >
class ConformalMatrixCoefficients : public MatrixCoefficients< TInputMesh >
{
public:
    typedef MatrixCoefficients< TInputMesh > Superclass;

    typedef TInputMesh InputMeshType;
    typedef typename InputMeshType::CoordRepType    InputCoordRepType;
    typedef typename InputMeshType::PointType       InputPointType;
    typedef typename InputMeshType::PointIdentifier InputPointIdentifier;
    typedef typename InputMeshType::QEType          InputQEType;

    ConformalMatrixCoefficients( ) :
        Superclass( )
    { };

    /**
    * \param[in] iMesh
    * \param[in] iEdge
    * \return \f$ \text{cot} \alpha_{ij} + \text{cot} \beta_{ij} \f$
     */
    InputCoordRepType operator ( ) ( InputMeshType* iMesh,
                                InputQEType* iEdge ) const
    {
        InputPointIdentifier id1 = iEdge->GetOrigin( );
        InputPointIdentifier id2 = iEdge->GetDestination( );
        InputPointType pt1 = iMesh->GetPoint( id1 );
        InputPointType pt2 = iMesh->GetPoint( id2 );
        
        InputCoordRepType oValue( 0. );

        if( iEdge->IsLeftSet() )
          {
          InputPointIdentifier idA = iEdge->GetLnext( )->GetDestination( );
          InputPointType ptA = iMesh->GetPoint( idA );
          oValue += Triangle< InputPointType >::Cotangent( pt1, ptA, pt2 );
          }
        if( iEdge->IsRightSet() )
          {
          InputPointIdentifier idB = iEdge->GetRnext( )->GetOrigin( );
          InputPointType ptB = iMesh->GetPoint( idB );
          oValue += Triangle< InputPointType >::Cotangent( pt1, ptB, pt2 );
          }

        return vnl_math_max( static_cast< InputCoordRepType >( 0. ), oValue );
    };
};

/**
    * \brief Compute a matrix filled with Authalic Coefiicients of the edge,
    *        wherever two vertices are connected with an edge.
    * \note  Belongs to the Parameterisation package.
    * \note  See paper:
    */
template< typename TInputMesh >
class AuthalicMatrixCoefficients : public MatrixCoefficients< TInputMesh >
{
public:
    typedef MatrixCoefficients< TInputMesh > Superclass;

    typedef TInputMesh InputMeshType;
    typedef typename InputMeshType::CoordRepType    InputCoordRepType;
    typedef typename InputMeshType::PointType       InputPointType;
    typedef typename InputMeshType::PointIdentifier InputPointIdentifier;
    typedef typename InputMeshType::QEType          InputQEType;

    AuthalicMatrixCoefficients( ) : Superclass( )
    { };

    /**
     * \param[in] iMesh
     * \param[in] iEdge
     * \return \f$ \frac{\text{cot} \gamma_{ij} + \text{cot}
     \delta_{ij}}{\|\boldsymbol{p1} - \boldsymbol{p2} \|} \f$
     */
    InputCoordRepType operator ( ) ( InputMeshType* iMesh,
                                    InputQEType* iEdge ) const
    {
        InputPointIdentifier id1 = iEdge->GetOrigin( );
        InputPointType pt1 = iMesh->GetPoint( id1 );
        
        InputPointIdentifier id2 = iEdge->GetDestination( );
        InputPointType pt2 = iMesh->GetPoint( id2 );
        
        InputCoordRepType oValue( 0. );

        if( iEdge->IsLeftSet() )
        {
          InputPointIdentifier idA = iEdge->GetLnext( )->GetDestination( );
          InputPointType ptA = iMesh->GetPoint( idA );
          oValue +=
            Triangle< InputPointType >::Cotangent( pt1, pt2, ptA );
        }
        if( iEdge->IsRightSet() )
        {
          InputPointIdentifier idB = iEdge->GetRnext( )->GetOrigin( );
          InputPointType ptB = iMesh->GetPoint( idB );
          oValue +=
            Triangle< InputPointType >::Cotangent( pt1, pt2, ptB );
        }

        return oValue / pt1.EuclideanDistanceTo( pt2 );
    };
};

/**
    * \brief Compute a mtrix filled by intrinsic Coefficients of the edge,
    *        wherever two vertices are connected by an edge.
    * \note  Belongs to the parameterization Package.
    * \note  See paper:
    */
template< typename TInputMesh >
class IntrinsicMatrixCoefficients : public MatrixCoefficients< TInputMesh >
{
public:
    typedef MatrixCoefficients< TInputMesh > Superclass;

    typedef TInputMesh InputMeshType;
    typedef typename InputMeshType::CoordRepType    InputCoordRepType;
    typedef typename InputMeshType::QEType          InputQEType;

    InputCoordRepType m_Lambda;

    IntrinsicMatrixCoefficients( const InputCoordRepType& iLambda ) :
        m_Lambda( iLambda )
    { };

    InputCoordRepType operator ( ) ( InputMeshType* iMesh,
                                    InputQEType* iEdge ) const
    {
        AuthalicMatrixCoefficients< TInputMesh > authalic;
        ConformalMatrixCoefficients< TInputMesh > conformal;

        InputCoordRepType oValue = m_Lambda * conformal( iMesh, iEdge )
            + ( 1. - m_Lambda ) * authalic( iMesh, iEdge );

        return oValue;
    };
};

/**
    * \brief Compute a matrix filled with Harmonic coefficients, wherever
    *        two vertices are connected by an edge.
    * \note  Belongs to the parameterization package.
    * \note  See paper:
    */
template< typename TInputMesh >
class HarmonicMatrixCoefficients : public MatrixCoefficients< TInputMesh >
{
public:
    typedef MatrixCoefficients< TInputMesh > Superclass;

    typedef TInputMesh InputMeshType;
    typedef typename InputMeshType::CoordRepType    InputCoordRepType;
    typedef typename InputMeshType::PointType       InputPointType;
    typedef typename InputPointType::VectorType InputVectorType;
    typedef typename InputMeshType::PointIdentifier InputPointIdentifier;
    typedef typename InputMeshType::QEType          InputQEType;

    itkStaticConstMacro( PointDimension, unsigned int,
        InputPointType::PointDimension );


    HarmonicMatrixCoefficients( ) :
        Superclass( )
    { };


    InputCoordRepType operator ( ) ( InputMeshType* iMesh,
                                InputQEType* iEdge ) const
    {
        InputPointIdentifier id1 = iEdge->GetOrigin( );
        InputPointIdentifier id2 = iEdge->GetDestination( );

        InputPointIdentifier idA = iEdge->GetLnext( )->GetDestination( );
        InputPointIdentifier idB = iEdge->GetRnext( )->GetOrigin( );

        InputPointType pt1 = iMesh->GetPoint( id1 );
        InputPointType pt2 = iMesh->GetPoint( id2 );
        InputPointType ptA = iMesh->GetPoint( idA );
        InputPointType ptB = iMesh->GetPoint( idB );

        InputVectorType v1A = ptA - pt1;
        InputVectorType v1B = ptB - pt1;
        InputVectorType v12 = pt2 - pt1;

        InputCoordRepType L1A = v1A * v1A;
        InputCoordRepType L1B = v1B * v1B;
        InputCoordRepType L12 = v12 * v12;

        InputCoordRepType L2A = pt2.SquaredEuclideanDistanceTo( ptA );
        InputCoordRepType L2B = pt2.SquaredEuclideanDistanceTo( ptB );

        Cross< InputVectorType > cross;

        InputCoordRepType AreaA = 0.5 * ( cross( v1A, v12 ).GetNorm( ) );
        InputCoordRepType AreaB = 0.5 * ( cross( v1B, v12 ).GetNorm( ) );

        InputCoordRepType
            oValue = ( L1A + L2A - L12 ) / AreaA + ( L1B + L2B - L12 ) / AreaB;

        return oValue;
    };
};

}
#endif
