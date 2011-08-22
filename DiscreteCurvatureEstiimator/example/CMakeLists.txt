INCLUDE_DIRECTORIES( BEFORE
  ${ITK_INCLUDE_DIR}
  "${PROJECT_SOURCE_DIR}/src/" )

SET( CURVATURE_EXAMPLE DiscreteCurvatureEstimator )

ADD_EXECUTABLE( ${CURVATURE_EXAMPLE} ${CURVATURE_EXAMPLE}.cxx )
TARGET_LINK_LIBRARIES( ${CURVATURE_EXAMPLE} ITKCommon ITKQuadEdgeMesh )   

FOREACH( curvature_type 0 1 2 3 )
  ADD_TEST( DiscreteCurvatureEstimatorTest${curvature_type}
      DiscreteCurvatureEstimator
      ${DATAROOT}/sphere.vtk ${curvature_type} )
ENDFOREACH( curvature_type )