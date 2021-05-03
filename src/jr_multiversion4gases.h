//   Since we cannot make use of template-supported 
//   multi-versioning kernel<int CO2, int H2O, int N2, int O2>
//   this workaround does the same in C-code.

#define  CAT(a,b) a##b
#define XCAT(a,b) CAT(a,b)

#define CO2 0
  #define H2O 0
    #define N2 0
      #define O2 0
          #include KERNEL
      #undef  O2
      
      #define O2 1
          #include KERNEL
      #undef  O2
    #undef  N2

    #define N2 1
      #define O2 0
          #include KERNEL
      #undef  O2

      #define O2 1
          #include KERNEL
      #undef  O2
    #undef  N2
  #undef  H2O

  #define H2O 1
    #define N2 0
      #define O2 0
          #include KERNEL
      #undef  O2

      #define O2 1
          #include KERNEL
      #undef  O2
    #undef  N2

    #define N2 1
      #define O2 0
          #include KERNEL
      #undef  O2

      #define O2 1
          #include KERNEL
      #undef  O2
    #undef  N2
  #undef  H2O
#undef  CO2

#define CO2 1
  #define H2O 0
    #define N2 0
      #define O2 0
          #include KERNEL
      #undef  O2

      #define O2 1
          #include KERNEL
      #undef  O2
    #undef  N2

    #define N2 1
      #define O2 0
          #include KERNEL
      #undef  O2

      #define O2 1
          #include KERNEL
      #undef  O2
    #undef  N2
  #undef  H2O
  
  #define H2O 1
    #define N2 0
      #define O2 0
          #include KERNEL
      #undef  O2

      #define O2 1
          #include KERNEL
      #undef  O2
    #undef  N2

    #define N2 1
      #define O2 0
          #include KERNEL
      #undef  O2

      #define O2 1
          #include KERNEL
      #undef  O2
    #undef  N2
  #undef  H2O
#undef  CO2

#undef XCAT
#undef  CAT
