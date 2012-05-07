void (*createRngArray)(int newNObj, int initialize);
void (*initializeRngInstance)(int k, unsigned int s);
void (*destroyRngArray)();
double (*getRngDouble)(int k);
void (*setRngParams)(int initType);

