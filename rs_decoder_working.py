import galois
import numpy as np
import time

verbose = True
GF = galois.GF(2**4)
if verbose: print(GF.repr_table())
alpha = GF.primitive_element


# ------------------------------------------------------------
# find_syndromes()
# 
# This function takes in two inputs:
# "received" is the polynomial that has been transmitted to us 
# "T" is the amount of errors we can correct
#
# The function can return two different outputs depending on the input:
# If syndromes are found, it returns a list of syndromes
# in the Galois field (S) and a list of syndromes as integers (S_powers)
#
# If all the syndromes are 0, we return -1 to signal that
# the codeword is correct
#
# It also makes use of the primitive element, 
# as defined globally as "alpha"
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# euclidean_algorithm()
#
# Takes in the polynomial that has been constructed
# using S_powers as returned from find_syndromes(), as well as
# the number of errorrs that can be corrected (T)
#
# After performing the Euclidean Algorithm, the function
# returns the quotient (q), the remainder (r) and the function (g)
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# error_correction()
#
# takes in 3 inputs:
# - rx is the remainder polynomial from euclidean_algorithm()
# - gx is the calculated polynomial from euclidean_algorithm()
# - T is the number of errors that can be corrected
#
# this function finds roots and error magnitude, and returns
# a polynomial that will correct any errors in the original message
# ------------------------------------------------------------
def error_correction(rx, gx, T):
    roots = gx.roots()
    int_roots = GF.log(roots)
    if verbose: print(f"Roots in GF: {roots}\nRoots in int: {int_roots}")

    gx_prime = gx.derivative()
    error_mag = []
    for i in range(len(roots)):
        error_mag.append(-roots[i]**(-(2*T+1)) * (rx(roots[i]) / gx_prime(roots[i])))

    error_vector = GF(np.zeros(15, dtype=int)) # this should be made dynamic
    error_vector[int_roots] = error_mag
    error_poly = galois.Poly(np.flip(error_vector), field=GF)

    if verbose: print(f"g'(x): {gx_prime}\nError magnitude: {error_mag}\nError vector: {error_vector}\nError polynomial: {error_poly}")

    return error_poly

# Correct codeword:
# Poly(x^8 + 7x^7 + 9x^6 + 3x^5 + 12x^4 + 10x^3 + 12x^2, GF(2^4))
#
# Codeword from the book:
# Poly(x^8 + 9x^6 + 3x^5 + 10x^3 + 12x^2 + 2)
# 
# Codeword as written in the book: (a = alpha)
# r(x) = x^8 + a^14 * x^6 + a^4 * x^5 + a^9 * x^3 + a^6 * x^2 + a

# CORRECT CODEWORD
#degrees = [8,7,6,5,4,3,2] 
#coeffs = [GF(1), GF(7), GF(9), GF(3), GF(12), GF(10), GF(12)]

# CODEWORD FROM THE BOOK (WITH 3 ERRORS)
degrees = [8, 6, 5, 3, 2, 0] 
coeffs = [GF(1), alpha**14, alpha**4, alpha**9, alpha**6, alpha**1]

# CODEWORD WITH RANDOM ERRORS
#degrees = [8,7,6,1,4,5,2] 
#coeffs = [GF(11), GF(5), GF(9), GF(3), GF(12), GF(10), GF(13)]

received = galois.Poly.Degrees(degrees, coeffs=coeffs, field=GF)
print(f"Received codeword: {received}")

T = 3

time_start = time.time()
S, S_powers = find_syndromes(received, T)
if S != -1:
    S_alpha = []
    for i in range(len(S_powers)):
        S_alpha.append(alpha**S_powers[i])
    
    S_poly = galois.Poly(S_alpha, field=GF)
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

else:
    print(f"Codeword already correct")

time_slut = time.time()
time_elapsed = time_slut - time_start
print(f"Decoding took {str(time_elapsed)[:4]} seconds")
input("Press enter to exit...")