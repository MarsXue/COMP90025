## COMP90025 Project 2
The focus of this project is to simulate the two-dimensional N-body problem.

### Team Member
- Hanbin Li - $hanbinl1$
- Wenqing Xue - $wenqingx$

### Implementation
- Basic OpenMP + MPI approach: $N^2$ complexity
  - N2_code.c
- Barnes-Hut OpenMP + MPI approach: $N\log N$ complexity
  - NlogN_code.c
  - tree.c

### Folder Structure
```js
/    
├── N2
│   └── N2_code.c       // C code for N^2 approach
├── NlogN
│   ├── NlogN_code.c    // C code for BH approach
│   └── tree.c          // C code for BH tree
├── README.md
├── figure.py           // Python code for figures
└── generator.c         // C code for generating the test data
```
