void (*createRngArray)(int new_n, int initialize);
void (*initializeRngInstance)(int k, unsigned int s);
void (*destroyRngArray)();
double (*getRngDouble)(int k);

