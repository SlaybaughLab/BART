#ifndef __in_group_base_h__
#define __in_group_base_h__

class InGroupBase : public IterationBase
{
public:
  InGroupBase ();
  virtual ~InGroupBase();
  
  // has to be overridden for NDA
  virtual void solve_in_group ();
};

#endif //__in_group_base_h__
