#ifndef __VTKPVORONOIREADER_h
#define __VTKPVORONOIREADER_h
 
#include <vtkMultiBlockDataSetAlgorithm.h>
#include <vtkSmartPointer.h>
#include "voronoi.h"

class vtkMultiProcessController;
class vtkUnstructuredGrid;

class vtkPVoronoiReader : public vtkMultiBlockDataSetAlgorithm
{
public:
  vtkTypeMacro(vtkPVoronoiReader,vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static vtkPVoronoiReader *New();
 
  // Description:
  // Specify file name of the .abc file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
 
protected:
  vtkPVoronoiReader();
  ~vtkPVoronoiReader(){}
 
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
 
private:
  vtkPVoronoiReader(const vtkPVoronoiReader&);  // Not implemented.
  void operator=(const vtkPVoronoiReader&);  // Not implemented.

  //ser_io
  void ReadFooter(FILE*& fd, int64_t*& ftr, int& tb);
  void ReadHeader(FILE *fd, int *hdr, int64_t ofst);
  int CopyHeader(unsigned char *in_buf, int *hdr);
  void ReadBlock(FILE *fd, vblock_t* &v, int64_t ofst);
  void CopyBlock(unsigned char *in_buf, vblock_t* &v);
  #if 0 // not using compression for now
    void DecompressBlock(unsigned  char* in_buf, int in_size, 
		       vector<unsigned char> *decomp_buf, 
		       int *decomp_size);
  #endif

  int dim; // number of dimensions in the dataset
  bool swap_bytes; // whether to swap bytes for endian conversion

  //swap
  void Swap(char *n, int nitems, int item_size);
  void Swap8(char *n);
  void Swap4(char *n);
  void Swap2(char *n);
 
  //parallelism
  int NumProcesses;
  int MyId;
  vtkMultiProcessController *Controller;
  void SetController(vtkMultiProcessController *c);

  //voronoi
  void vor2ugrid(struct vblock_t *block, vtkSmartPointer<vtkUnstructuredGrid> &ugrid);
 
  char* FileName;
};
 
#endif
