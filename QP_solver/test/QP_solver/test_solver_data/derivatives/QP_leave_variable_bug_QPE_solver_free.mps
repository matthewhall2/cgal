* Number-type: integer
* Description: Freed instance of original file
* Generated-by: master_mps_to_derivatives-create_free_instance
NAME QP_leave_variable_bug_QPE_solver
ROWS
  N obj
  G c0
  G c1
  G c2
  G c3
  G c4
COLUMNS
  x0  obj  0
  x0  c0  -4
  x0  c1  1
  x0  c2  1
  x0  c3  0
  x0  c4  0
  x1  obj  5
  x1  c0  2
  x1  c1  1
  x1  c2  0
  x1  c3  1
  x1  c4  0
  x2  obj  0
  x2  c0  0
  x2  c1  1
  x2  c2  0
  x2  c3  0
  x2  c4  1
RHS
  rhs c0  -8
  rhs c1  2
  rhs c2  0
  rhs c3  0
  rhs c4  0
BOUNDS
  MI  BND  x0
  MI  BND  x1
  MI  BND  x2
QMATRIX
  x0  x0  128
  x0  x1  -32
  x1  x0  -32
  x1  x1  8
  x1  x2  4
  x2  x1  -4
ENDATA
