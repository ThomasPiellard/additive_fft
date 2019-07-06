'''
Construction of constant arithmetic objects
'''

def next_g(prev_f, z, w, o_ring, base_ring):
    '''
        # brief
            constructing recursively the polynomials g and f

        # Args
            prev_f previous f in base field polynomial
            z generator of the polyring over smallest ext containing primitive 2P-1 roots fo unitty
            w 2p-1 root prim root of unity
    '''
    charac = o_ring.characteristic() 
    next_g = 1
    _prev_f = o_ring(prev_f)
    for i in range(0,2*charac-1):
        next_g *= _prev_f(z * (w**i))
    cont = next_g.dict()
    _next_g = 0
    for k in cont.keys():
        _next_g += cont[k]*(z**(k/(2*charac-1)))
    _next_g = base_ring(_next_g)
    return _next_g

def next_f(prev_g, x, charac):
    '''
        # brief
            generates the next f from preivous g (that is supposed to be irr)
            x generator of the polyring over the base field

        # args
            prev_g previous g (expressed in polyring over base field)
    '''
    next_f = prev_g(x**charac - x)
    return next_f

def compute_ith_f(f_1, x, z, o_ring, base_ring, w, nb_iter):
    '''
        # brief
            computes the i-th polynomial f_i. If i is a power of the charac, the polynomial is irreducible

        # params
            f_1 the initial polynomial x**p-x-1
            x the generator of the polynomial ring over the base field
            z the generator of the polynomial ring over the field containing w
            w the primitive root of 2*charac-1
            nb_iter the number of iterations to go through
    '''
    charac = base_ring.characteristic() 
    current_g = 0
    current_f = f_1
    res = [current_f]
    for i in range(0, nb_iter):
        current_g = next_g(current_f, z, w, o_ring, base_ring)
        current_f = next_f(current_g, x, charac)
        res += [current_f]
    return res

def setup_basis(charac, degree):
    '''
    brief
        setup the parameters for the additive fft.
        The basis (y_0,..,y_m) verifies q(y_m)=y_(m-1)

        Elements are encoded as base charac numbers: a_i*charac**i, a_i <charac, i<charac**degree
    params
        charac the characteristic of the base field
    degree degree of the extension (the degree of the extension is charac**degree)
    returns
        tuple (field, basis)
    ''' 
    base_field = GF(charac)
    poly_base_field.<x> = PolynomialRing(base_field) 

    # computation of w (primitive 2charac -1 root of 1) awa the smallest field in which it belongs
    smallest_extension=1
    while (charac**smallest_extension-1)%(2*charac-1) != 0:
       smallest_extension+=1 

    omega_field = GF(charac**smallest_extension, 'a', modulus="primitive")
    a = omega_field.gen()
    w = a**((charac**smallest_extension-1)/(2*charac-1))

    omega_ring.<z> = PolynomialRing(omega_field)
    z = omega_ring.gen()

    # computation of the sequence of polynomials generating intermediate fields
    f_1 = x**charac -x - 1
    f_s = compute_ith_f(f_1, x, z, omega_ring, poly_base_field, w, degree-1)
    big_field = GF(charac**(charac**degree), 'a', modulus = f_s[-1])
    
    # getting back the f_is in the biggest polynomial ring
    R.<u> = PolynomialRing(big_field) 

    for i,f in enumerate(f_s):
        d = f.dict()
        tmp = 0
        for k in d.keys():
            tmp += big_field(d[k])*u**k
        f_s[i] = tmp

    # compute the u_i
    u_i = [f.roots()[0][0] for f in f_s]
     
    # computation of the basis
    tmp = 1
    for i in u_i:
        tmp = tmp*i**(charac-1)
    
    basis = [tmp]
    comp = u**charac - u
    for i in range(charac**degree - 1):
        basis += [comp(basis[-1])]
    basis.reverse()

    return big_field, basis

def additive_fft(f, basis, field):
    '''
    brief
        computes the additive fft of f, on the space spanned by basis

        This version mimics exactly the multiplicative (with roots of unity) Fourier transform.
        The homomorphism here is x**charac - x (instead of x**2 in the multiplicative version)
        
        In one call, we compute Q(u,v) = P(u) mod (u**charac - u - v) (in lexicographic order in the bivariate ring)
        We end up with charac polynomials in v: f_0,..,f_v-1
        then: for v in (u**charac - u):
                res += sum u**i * fft(f_i) for i in 0..charac-1 for the charac u s.t. u**charac -u = v 

    params
        f function to evaluate
        basis the basis describing the vector space (y1_,..,y_m with S(y_j) = y_j-1)
        field the field containing the vector space
    returns 
        a list of containing the evaluation of f
    '''
    charac = field.characteristic()
    R.<x> = PolynomialRing(field)
    q = x**charac - x
    t = time.time()

    def rec_additive_fft(_f, basis):
       
        res = [0]*charac**len(basis) 
       
        if len(basis)==1: # this depends on q
            return [_f(i*basis[0]) for i in range(charac)] 
            
        else:
            extracted_poly = extract_poly(_f, q, x)
            _int_res = [rec_additive_fft(poly, basis[:-1]) for poly in extracted_poly]
            
            # for each element in Im(q) 
            for i in range(charac**(len(basis)-1)):
              
                tmp = i 
                b = field(0) 
                for j in range(len(basis)-1):
                    b += (tmp%charac)*basis[j+1]
                    tmp = tmp //charac
               
                tmp = charac*i 

                # for each x in preimage of i 
                for k in range(charac): # depends on q
                    
                    # translate into field element: 
                    _b = b + k*basis[0] 

                    # final computation using x
                    for m,n in enumerate(_int_res): 
                        res[tmp]+=n[i]*_b**m

                    tmp+=1
                    
            return res
    
    res = rec_additive_fft(f, basis)

    t = time.time() - t
    print('elapsed time: {}'.format(t))

    return res

def compute_sm(poly_ring, m):
    '''
    brief
        Computes the series of Sm up to m
    '''
    res = []
    x = poly_ring.gen()
    charac = poly_ring.characteristic()
    t =time.time()
    res += [x]

    for i in range(1,m): 
        tmp = 0
        for j in range(i+1):
            #tmp += (-1)**(i - j) * (binomial(i,j)%charac)*(x**(charac**j))
            tmp += (binomial(i,j)%charac)*x**(charac**j) 
        res += [tmp]
    t = time.time() - t
    print('elapsed time: {}'.format(t))
    return res

def compute_naive_sm(poly_ring, m):
    res = []
    x = poly_ring.gen()
    s = x**poly_ring.characteristic() - x
    t =time.time()
    res += [x]
    for i in range(m-1):
        res += [s(res[-1])]
    t = time.time() - t
    print('elapsed time: {}'.format(t))
    return res

def add_fft(f, basis, field):
    '''
    brief
        Other method for additive fft. Requires precomputations, but once those are done it is faster than the more
        generic one above.

    params
        f function to evaluate
        basis generator set of the vector space on which f is evaluated
        field the smallest field containing <basis>

    returns
        a list of evalution of f on <basis>, indexed as follows:
            i-th elemts is a0*basis[0]+..+an*basis[n] where ai is the decomposition of i in basis base charac.
    '''

    charac = field.characteristic()

    R.<x> = PolynomialRing(field) 
    q = x**charac - x
    tmp = q
    q_s = [x]

    #precomputation: takes A LOT of time!!! In non python, there is a way to speed things up big time
    for i in range(charac**len(basis) - 1):
        q_s += [q(q_s[-1])]
    print(q_s)
    
    #time is taken after the precomputations
    t = time.time()
        
    def rec_add_fft(f, mod, ind):
        
        if ind<0: 
            return [f] 

        else:
            el = 0
            for i,j in enumerate(mod):
                el += j*basis[i+1]
            res = []
            for i in range(charac):
                res += rec_add_fft(f % (q_s[ind] - (el+i)), [i]+mod, ind-1)
            return res

       
    r = rec_add_fft(f, [], len(basis)-1)
    t = time.time() - t
    print('elapsed time: {}'.format(t))
    return r

''' utils functions '''

def extract_poly(f, s, x):
    '''
    brief
        compute f % (s - v)
    '''
    t = [0*x] * s.degree()
    for i in range(f.degree()//s.degree() + 1):
        r = f % s
        d = r.dict()
        for k in d.keys():
            t[k] += d[k]*x**i
        f = f//s
    return t

def total_elts(basis, charac):
    '''
    brief
        compute all elements of the vector space generated by basis over a finite field of charac charac
    '''
    res = []
    for i in range(charac**len(basis)):
        tmp = i
        b = 0
        for k in range(len(basis)):
            b += (tmp%charac)*basis[k]
            tmp = (tmp - tmp%charac)/charac
        res += [b]
    return res

def naive_evaluation(f, basis, charac):
    '''
    brief
        computes the naive evaluation of f on the the vector space generated by basis over a finite
        field of characteristic charac
    '''
    elts = total_elts(basis, charac)
    
    t = time.time()
    r = []
    for i in elts:
        r += [f(i)]
    t = time.time() - t
    print('elapsed time: {}'.format(t))
    return r

def _binomial(a, b, m, charac):
    '''
    brief
        trick to compute binomial(a, b) in characteristic not 0
        Uses Lucas congruences (cf wiki)
    params
        a,b paramaters to compute (binomial(a, b))
        m power of charac bounding abov a & b
        charac characteristic
    '''
    if b>a:
        return 0
    elif b==a or b==0:
        return 1
    else:
        tmp = 1
        for i in range(m):
            tmp *= binomial(a%charac, b%charac)%3
            if tmp==0 or b==0:
                break
            a = (a - a%charac)/charac
            b = (b - b%charac)/charac
        return tmp%charac

