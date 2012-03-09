class ClassRandom {
private:
    // internal state of generator mrg32k5a from L'Ecuyer99
    double s10, s11, s12, s13, s14, s20, s21, s22, s23, s24;
    double MRG32k5a();

public:
    // constructors
    ClassRandom() {
        s10 = s11 = s12 = s13 = s14 = 12345.0;
        s20 = s21 = s22 = s23 = s24 = 12345.0;
    }
    ClassRandom(double *seed) {
        putState(seed);
    }
    ~ClassRandom() {};
    // interface
    inline void putState(double *state) {
        s10 = state[0];
        s11 = state[1];
        s12 = state[2];
        s13 = state[3];
        s14 = state[4];
        s20 = state[5];
        s21 = state[6];
        s22 = state[7];
        s23 = state[8];
        s24 = state[9];
    }
    inline double getDouble() {
        return MRG32k5a();
    }
};

typedef ClassRandom *PtrRandom;

