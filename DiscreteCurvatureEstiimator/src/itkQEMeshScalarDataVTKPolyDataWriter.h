#ifndef __itkQEMeshScalarDataVTKPolyDataWriter_h
#define __itkQEMeshScalarDataVTKPolyDataWriter_h

#include <itkVTKPolyDataWriter.h>
#include <fstream>

namespace itk
{
  /**
   * \class QEMeshScalarDataVTKPolyData
   * \brief
  */
  template< class TMesh >
  class QEMeshScalarDataVTKPolyDataWriter : public VTKPolyDataWriter< TMesh >
  {
  public:
    typedef QEMeshScalarDataVTKPolyDataWriter Self;
    typedef SmartPointer< Self > Pointer;
    typedef SmartPointer< const Self > ConstPointer;
    typedef VTKPolyDataWriter< TMesh > Superclass;

    /** Run-time type information (and related methods).   */
    itkTypeMacro( QEMeshScalarDataVTKPolyDataWriter, VTKPolyDataWriter );

    /** New macro for creation of through a Smart Pointer   */
    itkNewMacro( Self );

    typedef TMesh MeshType;
    typedef typename MeshType::Pointer MeshPointer;
    typedef typename MeshType::CellType CellType;

    typedef typename MeshType::PointsContainerPointer PointsContainerPointer;
    typedef typename MeshType::PointsContainerIterator PointsContainerIterator;

    typedef typename MeshType::PointDataContainerPointer
      PointDataContainerPointer;
    typedef typename MeshType::PointDataContainerIterator
      PointDataContainerIterator;

    typedef typename MeshType::CellsContainerPointer CellsContainerPointer;
    typedef typename MeshType::CellsContainerIterator CellsContainerIterator;

    typedef typename MeshType::CellDataContainerPointer CellDataContainerPointer;
    typedef typename MeshType::CellDataContainerIterator
      CellDataContainerIterator;
    
    /** Set/Get the name of the CellDataName where data are written. */
    itkSetStringMacro(CellDataName);
    itkGetStringMacro(CellDataName);

    /** Set/Get the name of the PointDataName where data are written. */
    itkSetStringMacro(PointDataName);
    itkGetStringMacro(PointDataName);

  protected:
    QEMeshScalarDataVTKPolyDataWriter() : Superclass(), m_CellDataName( "" ),
      m_PointDataName( "" ) {}
    ~QEMeshScalarDataVTKPolyDataWriter() {}

    std::string m_CellDataName;
    std::string m_PointDataName;

    void GenerateData()
    {
      Superclass::GenerateData();
      WriteCellData();
      WritePointData();
    }

    void WriteCellData()
    {
      CellDataContainerPointer celldata = this->m_Input->GetCellData();

      if( celldata )
        {
        if( celldata->Size() != 0 )
          {
          std::ofstream outputFile( this->m_FileName.c_str(),
            std::ios_base::app );
          outputFile <<"CELL_DATA " <<this->m_Input->GetNumberOfFaces()
            <<std::endl;
          outputFile <<"SCALARS ";
          if( m_CellDataName != "" )
            outputFile <<m_CellDataName <<" " <<m_CellDataName <<std::endl;
          else
            outputFile <<"double double" <<std::endl;
          
          outputFile <<"LOOKUP_TABLE default" <<std::endl;

          unsigned long k(0);

          CellsContainerPointer cells = this->m_Input->GetCells();
          CellsContainerIterator it = cells->Begin();

          for( CellDataContainerIterator c_it = celldata->Begin();
              c_it != celldata->End();
              ++c_it, ++it )
            {
            CellType* cellPointer = it.Value();
            if( cellPointer->GetType() != 1 )
              {
              outputFile <<c_it.Value();
              if( k++ % 3 == 0 )
                outputFile <<std::endl;
              }
            }
          outputFile <<std::endl;
          outputFile.close();
          }
        }
    }

    void WritePointData()
    {
      PointDataContainerPointer pointdata = this->m_Input->GetPointData();

      if( pointdata )
        {
        std::ofstream outputFile( this->m_FileName.c_str(), std::ios_base::app );
        outputFile <<"POINT_DATA " <<this->m_Input->GetNumberOfPoints()
          <<std::endl;
        outputFile <<"SCALARS ";
        if( m_PointDataName != "" )
          outputFile <<m_PointDataName <<" " <<m_PointDataName <<std::endl;
        else
          outputFile <<"double double"<<std::endl;

        outputFile <<"LOOKUP_TABLE default" <<std::endl; 
        unsigned long k(0);

        for( PointDataContainerIterator c_it = pointdata->Begin();
            c_it != pointdata->End();
            ++c_it, ++k )
          {
          outputFile <<c_it.Value() <<" ";
          if( k % 3 == 0 )
            outputFile <<std::endl;
          }
        outputFile <<std::endl;
        outputFile.close();
        }
    }

  private:
    QEMeshScalarDataVTKPolyDataWriter( const Self& );
    void operator = ( const Self& );
  };
}
#endif
