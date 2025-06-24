import galois
import numpy as np

GF = galois.GF(2**4)
print(GF.repr_table())
alpha = GF.primitive_element
verbose = True

def find_syndromes(received, T):
    S = []
    for i in range(1, 2*T+1):
        Si = received(alpha**i)
        S.append(Si)
    if np.count_nonzero(S) == 0:
        return -1, -1
    else:
        S_powers = []
        for s in S:
            S_powers.append(GF.log(s))
        if verbose: print(f"Syndromes: {S_powers}")
        return S, S_powers

def euclidean_algorithm(S_poly, T):
    ax = galois.Poly.Degrees([2*T],field=GF)
    bx = S_poly
    r_prev = ax
    r = r_prev
    g_prev = GF(0)

    r_0 = bx
    g_0 = GF(1)

    if verbose: print(f"""Euclidean algorithm:\n-------------\nInitial variables:\n-------------\nq-1: -\nr-1: {r_prev}\ng-1: {g_prev}\n------------\nq0: -\nr0: {r_0}\ng0: {g_0}\n-------------""")
    i = 1
    while r.degree > T-1:
        q, r = divmod(r_prev, r_0)
        g = g_prev - (q * g_0)

        if verbose: print(f"""q{i}: {q}\nr{i}: {r}\ng{i}: {g}\n-------------""")
        g_prev = g_0
        g_0 = g
        r_prev = r_0
        r_0 = r
        i += 1
    
    return q, r, g

def error_correction(rx, gx, T):
    roots = gx.roots()
    int_roots = GF.log(roots)
    if verbose: print(f"Roots in GF: {roots}\nRoots in int: {int_roots}")

    gx_prime = gx.derivative()
    error_mag = [] # error magnitude
    for i in range(len(roots)):
        error_mag.append(-roots[i]**(-(2*T+1)) * (rx(roots[i]) / gx_prime(roots[i])))

    error_vector = GF(np.zeros(15, dtype=int)) ######## NEEDS TO BE MADE DYNAMIC
    error_vector[int_roots] = error_mag
    error_poly = galois.Poly(np.flip(error_vector), field=GF)

    if verbose: print(f"g'(x): {gx_prime}\nError magnitude: {error_mag}\nError vector: {error_vector}\nError polynomial: {error_poly}")

    return error_poly

# Poly(x^8 + 7x^7 + 9x^6 + 3x^5 + 12x^4 + 10x^3 + 12x^2, GF(2^4)) KORREKTE KODEORD

#degrees = [8,7,6,5,4,3,2] # KORREKTE KODEORD
#coeffs = [GF(1), GF(7), GF(9), GF(3), GF(12), GF(10), GF(12)]

#degrees = [8,7,6,1,4,5,2] # MED FEJL
#coeffs = [GF(11), GF(5), GF(9), GF(3), GF(12), GF(10), GF(13)]

degrees = [8, 6, 5, 3, 2, 0] # MED 3 FEJL (FRA BOG)
coeffs = [GF(1), alpha**14, alpha**4, alpha**9, alpha**6, alpha**1]

received = galois.Poly.Degrees(degrees, coeffs=coeffs, field=GF)
print(f"Received codeword: {received}")

T = 3

S, S_powers = find_syndromes(received, T)
if S == -1:
    print(f"Codeword already correct")
else:
    s_alpha = [alpha**S_powers[0],alpha**S_powers[1],alpha**S_powers[2],alpha**S_powers[3],alpha**S_powers[4],alpha**S_powers[5]]
    S_poly = galois.Poly(s_alpha, field=GF)
    if verbose: print(f"Syndrome polynomial (S(x)): {S_poly}")

    q, r, g = euclidean_algorithm(S_poly, T)

    error_poly = error_correction(r, g, T)

    corrected_codeword = error_poly + received
    print(f"Corrected codeword: {corrected_codeword}")

    S_correct = []
    for i in range(1, 2*T+1):
        Si = corrected_codeword(alpha**i)
        S_correct.append(Si)
    print(f"Syndromes for corrected codeword: {S_correct}")

    if np.count_nonzero(S_correct) == 0:
        print(f"Codeword has been decoded successfully")
    else:
        print(f"The codeword has not been correctly decoded")