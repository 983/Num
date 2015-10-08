#include "num.hpp"
#include "rat.hpp"
#include <assert.h>
#include <stdio.h>

int cmp(int a, int b){
    if (a == b) return 0;
    return a < b ? -1 : +1;
}

int cmp_abs(int a, int b){
    return cmp(std::abs(a), std::abs(b));
}

int gcd(int a, int b){
    a = std::abs(a);
    b = std::abs(b);
    while (1){
        if (b == 0) return a;
        a %= b;
        if (a == 0) return b;
        b %= a;
    }
}

void calculate_pqt(const Num &a, const Num &b, Num &P, Num &Q, Num &T){
    if (a + 1 == b){
        if (a.size()){
            P = (a*2 - 1)*(a*6 - 5)*(a*6 - 1);
            Q = a*a*a*"10939058860032000";
        }else{
            P = 1;
            Q = 1;
        }
        T = P * (a*545140134 + 13591409);
        if (a.get_bit(0)) T = -T;
    }else{
        Num P1, Q1, T1, P2, Q2, T2, m = (a + b) >> 1;
        calculate_pqt(a, m, P1, Q1, T1);
        calculate_pqt(m, b, P2, Q2, T2);
        P = P1 * P2;
        Q = Q1 * Q2;
        T = Q2 * T1 + P1 * T2;
    }
}

Num calculate_pi(size_t n_digits){
    double digits_per_iteration = 14.1816474627;
    Num iterations = n_digits / digits_per_iteration + 1;
    iterations += 3;
    Num P, Q, T;
    calculate_pqt(0, iterations, P, Q, T);
    Num pi = (Num(100).pow(n_digits)*10005).sqrt()*Q*426880 / T;
    return pi;
}

Rat calculate_e(int iterations = 10){
    Num x = 2 + 4*iterations;
    Rat result = x + 4;
    for (; x >= 6; x -= 4) result = result.inv() + x;
    return (result.inv() + 1).inv()*2 + 1;
}

int main(){
    Num pi = calculate_pi(1000);
    assert(pi == "31415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989");

    std::vector<char> c;
    pi.print(c);
    printf("pi to %i digits:\n%c.%s\n", int(c.size()) - 2, c[0], &c[1]);

    assert(pi == Num(&c[0]));

    Rat e = calculate_e();
    printf("\n15 digits of e: %.15f\n", e.to_double());
    double d = e.to_double() - 2.718281828459045;
    double eps = 1e-15;
    assert(-eps < d && d < eps);

    Rat r(1337, 191);
    r.print(c);
    assert(r == 7);

    assert(Num(666).mod_pow(1337, 1234) == 816);

    assert(Num(666).mod_pow(0, 1000) == 1);

    int n = 100;
    for (int a = -n; a <= n; a++){
        for (int b = -n; b <= n; b++){
            assert(cmp_abs(a, b) == Num::cmp_abs(Num(a), Num(b)));
            assert(cmp    (a, b) == Num::cmp(Num(a), Num(b)));
            assert((Num(a) == Num(b)) == (a == b));
            assert((Num(a) != Num(b)) == (a != b));
            assert((Num(a) <= Num(b)) == (a <= b));
            assert((Num(a) >= Num(b)) == (a >= b));
            assert((Num(a) <  Num(b)) == (a <  b));
            assert((Num(a) >  Num(b)) == (a >  b));
            assert((Num(a) +  Num(b)) == (a +  b));
            assert((Num(a) -  Num(b)) == (a -  b));
            assert((Num(a) *  Num(b)) == (a *  b));
            if (a >= 0 && b >= 0){
                assert(Num(a) >> b % 18 == Num(a >> b % 18));
                assert(Num(a) << b % 18 == Num(a << b % 18));
            }
            if (b != 0){
                assert((Num(a) /  Num(b)) == (a /  b));
                assert((Num(a) %  Num(b)) == (a %  b));
                if (a != 0) assert(Num::gcd(a, b) == gcd(a, b));
            }
            if (b > 1 && std::abs(b) <= Num::word_half_mask()){
                Num tmp(a);
                Num::word remainder;
                Num::div_mod_half_word(tmp, b, tmp, remainder);
                assert(remainder == std::abs(a%b));
                assert(tmp == Num(a/b));
            }
        }
    }

    for (size_t n = 2; n < 100; n++){
        Num a(n, Num::word_mask());
        Num c = a * a;
        size_t m = c.size();
        assert(c[0] == 1);
        for (size_t i = 1; i < m/2; i++) assert(c[i] == 0);
        assert(c[m/2] == Num::word_mask() - 1);
        for (size_t i = m/2 + 1; i < m; i++) assert(c[i] == Num::word_mask());
    }

	return 0;
}
