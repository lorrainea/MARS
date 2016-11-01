#define INFINITE 10000000.0 
#define NEG_INFTY -10000000.0

typedef struct Path_ {
        int *Minimum ;
        int *Maximum ; } Path ;

typedef struct Limits_ {
        Path Left ;
        Path Right ; } Limits ;

extern float RestrictedLevenshtein(int FirstCharacter, unsigned char *pattern1, unsigned char *pattern2, int length1, int length2, Limits Limit , Path *path, float *gamma, int *SlopeToMinimumPath) ;
void AllocateMemoryPath(Path *path, int size);
void FreeMemoryPath(Path *path);
extern void LimitInitialize(Limits *Limit, int length1, int length2) ;
extern float bb(float distance2,int length1, int length2, unsigned char *pattern1, unsigned char *pattern2, Path *BestPath, Path *path, int bb_type, float *gamma, int *cyc, unsigned int * rotation, unsigned int * distance, int *SlopeToMinimumPath) ;
extern void CopyPath(Path *Original, Path *Copy, int length) ;
float BoundFunction(int left, int right, float left_cost, float right_cost, unsigned char *pattern1, unsigned char *pattern2, int length1, int length2, int bound_type, float *gamma) ;
