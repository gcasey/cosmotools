#ifndef MERGERTREEFILEFORMAT_H_
#define MERGERTREEFILEFORMAT_H_


namespace cosmotk {

class MergerTreeFileFormat
{
public:
  enum
    {
    DIY = 0,
    GENERIC_IO_MPI,
    GENERIC_IO_POSIX,

    UNDEFINED
    };
};

}

#endif /* MERGERTREEFILEFORMAT_H_ */
