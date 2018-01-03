from scipy.misc import imread
from matplotlib import pyplot as plt
import pywt
import numpy as np
import scipy as sp
import Queue
import bitstring as bs

class WSQ:
    '''
    Perform compression using the Wavelet Scalar Quantization algorithm.
    All class attributes are set to None in __init__, but their values
    are initialized in the compress method.
    
    Attributes
    ----------
    _pixels : int, number of pixels in source image
    _s : float, scale parameter for image preprocessing
    _m : float, shift parameter for image preprocessing
    _Q : numpy array, quantization parameters q for each subband
    _Z : numpy array, quantization parameters z for each subband
    _bitstrings : list of 3 BitArrays, giving bit encodings for each group
    _tvals : tuple of 3 lists of bools, indicating which subbands in each
    groups were encoded
    _shapes : tuple of 3 lists of tuples, giving shapes of each subband in
    each group
    _huff_maps : list of 3 dicts, mapping huffman index to bit pattern
    '''
    
    def __init__(self):
        self._pixels = None
        self._s = None
        self._m = None
        self._Q = None
        self._Z = None
        self._bitstrings = None
        self._tvals = None
        self._shapes= None
        self._huff_maps = None
        self._img = None
        
    def _preProcess(self):
        m = self._img.mean()
        s = max(self._img.max()-m, m-self._img.min())/128
        self._m = m
        self._s = s
        self._img = (self._img - self._m )/ self._s

    def _postProcess(self):
        return self._img * self._s + self._m
        
    def _decompose16(self, img, wavelet):
        '''
        Decompose an image into 16 subbands
        '''
        subbands = []
        LL, HVD = pywt.dwt2(img, wavelet, mode='per')
        dec = pywt.dwt2(LL, wavelet, mode='per')
        subbands.append(dec[0])
        subbands.extend(dec[1])
        for i in xrange(3):
            dec = pywt.dwt2(HVD[i], wavelet, mode='per')
            subbands.append(dec[0])
            subbands.extend(dec[1])
        return subbands
    
    def _recreate16(self, subbands, wavelet):
        '''
        Recreate the image from the 16 subbands
        '''
        LL = pywt.idwt2((subbands[0], tuple(subbands[1:4])), wavelet, mode='per')
        details = []
        for i in xrange(1,4):
            details.append(pywt.idwt2((subbands[4*i], tuple(subbands[4*i+1:4*i+4])), wavelet, mode='per'))
        return pywt.idwt2((LL, tuple(details)), wavelet, mode='per')
    
    def _decompose(self):
        # Choose the type of wavelet to use
        wavelet = 'coif1'
        # Initialize subbands list so we have something to append to
        subbands = []
        # Decompose into initial 16 subbands
        temp1 = self._decompose16(self._img, wavelet)
        
        # Steps to decompose top left subband and its bottom and left 
        # neighbors into sixteen subbands apiece
        temp2 = []
        for i in xrange(3):
            temp2.append(self._decompose16(temp1[i], wavelet))
            
        # Decompose top left corner into four subbands
        ll, hvd = pywt.dwt2(temp2[0][0], wavelet, mode='per')
        
        # Append top left corner, 0
        subbands.append(ll)
        # Append the rest of the top left corner, 1-3
        subbands.extend(hvd)
        # Append the remaining part of that block (ignore what we
        # have already appended)
        subbands.extend(temp2[0][1:])
        # Append its right neighbor
        subbands.extend(temp2[1])
        # Append its bottom neighbor
        subbands.extend(temp2[2])
        # Append the rest of the original sixteen subbands
        subbands.extend(temp1[3:])
        return subbands
    
    def _recreate(self, subbands):
        # Choose the wavelet to use
        wavelet = 'coif1'
        # Put the tiny top left corner back together
        ll = pywt.idwt2((subbands[0], tuple(subbands[1:4])), wavelet, mode='per')
        
        temp1 = []
        temp2 = []
        
        # Put together tiny the top left subbblock
        temp2.append(ll)
        temp2.extend(subbands[4:19])
        # Merge together the new top left subblock
        temp1.append(self._recreate16(temp2, wavelet))
        # Merge its right neighbor
        temp1.append(self._recreate16(subbands[19:35], wavelet))
        # Merge its bottom neighbor
        temp1.append(self._recreate16(subbands[35:51], wavelet))
        # Add the rest of the subbands. Now we have the image decomposed
        # into sixteen subbands
        temp1.extend(subbands[51:])
        
        # Finish putting the image back together
        self._img = self._recreate16(temp1, wavelet)
        
    def _getBins(self, subbands, r, gamma):
        '''
        Calculate quantization bin widths for each subband.
        '''
        subband_vars = np.zeros(64)
        fracs = np.zeros(64)
        for i in xrange(len(subbands)):
            X, Y = subbands[i].shape
            fracs[i] = (X * Y) / (np.float(finger.shape[0]*finger.shape[1]))
            x = np.floor(X / 8.0)
            y = np.floor(9 * Y /32.0)
            Xp = np.floor(3 * X / 4.0)
            Yp = np.floor(7 * Y / 16.0)
            mu = subbands[i].mean()
            sigsq = (Xp * Yp -1)**(-1) * ((subbands[i][x:x+Xp, y:y+Yp] - mu)**2).sum()
            subband_vars[i] = sigsq
            
        A = np.ones(64)
        A[52], A[56] = [1.32]*2
        A[53], A[58], A[55], A[59] = [1.08]*4
        A[54], A[57] = [1.42]*2
        
        Qprime = np.zeros(64)
        mask = subband_vars >= 1.01
        Qprime[mask] = 10.0 / (A[mask]*np.log(subband_vars[mask]))
        Qprime[:4] = 1
        Qprime[60:] = 0
        
        K = []
        for i in xrange(60):
            if subband_vars[i] >= 1.01:
                K.append(i)
        
        while True:
            S = fracs[K].sum()
            P = ((np.sqrt(subband_vars[K]) / Qprime[K])**fracs[K]).prod()
            q = (gamma**(-1)) * (2**(r/S - 1)) * (P**(-1.0/S))
            E = []
            for i in K:
                if Qprime[i]/q >= 2 * gamma * np.sqrt(subband_vars[i]):
                    E.append(i)
            if len(E) > 0:
                for i in E:
                    K.remove(i)
                continue
            break
        
        Q = np.zeros(64)
        for i in K:
            Q[i] = Qprime[i]/q
        Z = 1.2*Q
        
        return Q, Z
        
    def _quantize(self, coeffs, Q, Z):
        '''
        Implement a uniform quantizer.
        
        Parameters
        ----------
        coeffs: numpy array containing the values to be quantized.
        Q: the step size of the quantization, a nonnegative float
        Z: the null-zone width (of the center/0 quantization bin)
                nonnegative float
        Returns
        -------
        numpy array of same shape as coeffs holding the quantized values
        '''
        if Q != 0:
            p = np.copy(coeffs)
            
            case1 = coeffs > Z/2.0
            case2 = np.abs(coeffs) <= Z/2.0
            case3 = coeffs < -Z/2.0
        
            p[case1] = np.floor((coeffs[case1]-Z/2.0)/Q) + 1
            p[case2] = 0
            p[case3] = np.ceil((coeffs[case3]+Z/2.0)/Q) - 1
        
        else:
            p = np.zeros(coeffs.shape)
        
        return p
    
    def _dequantize(self, coeffs, Q, Z, C=0.44):
        '''
        Approximately reverse the quantization
        
        Parameters
        ----------
        coeffs: numpy array of integer values
        Q, Z: see documentation for _quantize()
        
        Returns
        -------
        dequantized coefficients (of same shape as coeffs)
        '''
        a = np.copy(coeffs)
        case1 = coeffs > 0
        case2 = coeffs == 0
        case3 = coeffs < 0
        
        a[case1] = (coeffs[case1] - C) * Q + Z/2.0
        a[case2] = 0
        a[case3] = (coeffs[case3] + C) * Q - Z/2.0

        return a
        
    def _group(self, subbands):
        '''
        Split the quantized subbands into 3 groups
        
        Parameters
        ----------
        subbands - list of 64 numpy arrays containing quantized coefficients
        
        Returns
        -------
        (g1, g2, g3) - tuple of lists of integers corresponding to groups 1, 2, and 3
        (s1, s2, s3) - tuple of lists of tuples giving the shapes of the subbands
        
        '''
        g1 = []
        s1 = []
        t1 = []
        for i in xrange(19):
            s1.append(subbands[i].shape)
            if subbands[i].any():
                g1.extend(subbands[i].ravel())
                t1.append(True)
            else:
                t1.append(False)
            
        g2 = []
        s2 = []
        t2 = []
        for i in xrange(19,52):
            s2.append(subbands[i].shape)
            if subbands[i].any():
                g2.extend(subbands[i].ravel())
                t2.append(True)
            else:
                t2.append(False)
            
        g3 = []
        s3 = []
        t3 = []
        for i in xrange(52,64):
            s3.append(subbands[i].shape)
            if subbands[i].any():
                g3.extend(subbands[i].ravel())
                t3.append(True)
            else:
                t3.append(False)
                
        return (g1,g2,g3), (s1,s2,s3), (t1,t2,t3)
    
    def _ungroup(self, gs, ss, ts):
        '''
        Re-create the subbands list structure from the three groups.
        
        Parameters
        ----------
        groups: tuple of form (g1, g2, g3)
        shapes: tuple of form (s1, s2, s3)
        See the docstring for _group
        
        Returns
        -------
        subbands: list of 64 numpy arrays
        '''
        g1, g2, g3 = gs[0], gs[1], gs[2]
        s1, s2, s3 = ss[0], ss[1], ss[2]
        t1, t2, t3 = ts[0], ts[1], ts[2]
        
        subbands1 = []
        i = 0
        for j, shape in enumerate(s1):
            if t1[j]:
                l = shape[0] * shape[1]
                subbands1.append(np.array(g1[i:i+l]).reshape(shape))
                i += l
            else:
                subbands1.append(np.zeros(shape))
        
        subbands2 = []
        i = 0
        for j, shape in enumerate(s2):
            if t2[j]:
                l = shape[0] * shape[1]
                subbands2.append(np.array(g2[i:i+l]).reshape(shape))
                i += l
            else:
                subbands2.append(np.zeros(shape))
            
        subbands3 = []
        i = 0
        for j, shape in enumerate(s3):
            if t3[j]:
                l = shape[0] * shape[1]
                subbands3.append(np.array(g3[i:i+l]).reshape(shape))
                i += l
            else:
                subbands3.append(np.zeros(shape))

        return subbands1 + subbands2 + subbands3
            
    def _huffmanIndices(self, coeffs):
        '''
        
        '''
        N = len(coeffs)
        i = 0
        inds = []
        extra = []
        freqs = np.zeros(254)
        
        # Sweep through the quantized coefficients
        while i < N:
            # First handles runs of zeros
            zero_count = 0
            while coeffs[i] == 0:
                zero_count += 1
                i += 1
                if i >= N:
                    break
            if zero_count > 0 and zero_count < 101:
                inds.append(zero_count - 1)
                freqs[zero_count - 1] += 1
            elif zero_count >= 101 and zero_count < 256:
                inds.append(104)
                freqs[104] += 1
                extra.append(zero_count)
            elif zero_count >= 256:
                inds.append(105)
                freqs[105] += 1
                extra.append(zero_count)
            if i >= N:
                break
            
            # Now handle nonzero coefficients
            if coeffs[i] > 74 and coeffs[i] < 256:
                inds.append(100)
                freqs[100] += 1
                extra.append(coeffs[i])
            elif coeffs[i] >= 256:
                inds.append(102)
                freqs[102] += 1
                extra.append(coeffs[i])
            elif coeffs[i] < -73 and coeffs[i] > -256:
                inds.append(103)
                freqs[103] += 1
                extra.append(abs(coeffs[i]))
            else:
                inds.append(179 + coeffs[i])
                freqs[179 + coeffs[i]] += 1
            i += 1
        
        return inds, freqs, extra

    def _indicesToCoeffs(self, indices, extra):
        coeffs = []
        j = 0
        for s in indices:
            if s < 100:
                coeffs.extend(np.zeros(s+1))
            elif s == 104 or s == 105:
                coeffs.extend(np.zeros(extra[j]))
                j += 1
            elif s in [100, 101, 102, 103]:
                coeffs.append(extra[j])
                j += 1
            else:
                coeffs.append(s-179)
        return coeffs
    
    class huffmanLeaf():
        def __init__(self, symbol):
            self.symbol = symbol
        def makeMap(self, huff_map, path):
            huff_map[self.symbol] = path
            
    class huffmanNode():
        def __init__(self, left, right):
            self.left = left
            self.right = right
        def makeMap(self, huff_map, path):
            self.left.makeMap(huff_map, path + '0')
            self.right.makeMap(huff_map, path + '1')
            
    def huffman(self, freqs):
        q = Queue.PriorityQueue()
        for i in xrange(len(freqs)):
            leaf = self.huffmanLeaf(i)
            q.put((freqs[i], leaf))
        while q.qsize() > 1:
            l1 = q.get()
            l2 = q.get()
            weight = l1[0] + l2[0]
            node = self.huffmanNode(l1[1], l2[1])
            q.put((weight, node))
        huff_map = dict()
        root = q.get()[1]
        root.makeMap(huff_map, '')
        return huff_map
    
    def _encode(self, indices, extra, huff_map):
        bits = bs.BitArray()
        j = 0
        for s in indices:
            bits.append('0b' + huff_map[s])
            if s in [104, 100, 101]:
                bits.append('uint:8={}'.format(int(extra[j])))
                j += 1
            elif s in [102, 103, 105]:
                bits.append('uint:16={}'.format(int(extra[j])))
                j += 1
        return bits
    
    def _decode(self, bits, huff_map):
        indices = []
        extra = []
        
        dec_map = {v:k for k, v in huff_map.items()}
        bits = bs.ConstBitStream(bits)
        
        i = 0
        pattern = ''
        while i < bits.length:
            pattern += bits.read('bin:1')
            i += 1
            
            if dec_map.has_key(pattern):
                indices.append(dec_map[pattern])
                
                if dec_map[pattern] in (100, 104):
                    extra.append(bits.read('uint:8'))
                    i += 8
                
                elif dec_map[pattern] == 101:
                    extra.append(-bits.read('uint:8'))
                    i += 8
                    
                elif dec_map[pattern] in (102, 105):
                    extra.append(bits.read('uint:16'))
                    i += 16
                
                elif dec_map[pattern] == 103:
                    extra.append(-bits.read('uint:16'))
                pattern = ''
                
        return indices, extra
    
    def compress(self, img, r, gamma=2.5):
        '''
        The main compression routine. It computes and stores the bitstring
        representation of the compressed image, along with other values needed
        for decompressions.
        
        Parameters
        ----------
        img: numpy array containing 8-bit integer pixel values
        r: float, the closer to zero, the higher compression ratio
        gamma: float, a parameter used in quantization
        '''
        self._img = img
        self._pixels = self._img.shape[0] * self._img.shape[1]
        self._preProcess()
        subbands = self._decompose()
        self._Q, self._Z = self._getBins(subbands, r, gamma)
        
        q_subbands = [self._quantize(subbands[i], self._Q[i], self._Z[i]) for i in xrange(64)]
        
        groups, self._shapes, self._tvals = self._group(q_subbands)
        
        huff_maps = []
        bitstrings = []
        for i in xrange(3):
            inds, freqs, extra = self._huffmanIndices(groups[i])
            huff_map = self.huffman(freqs)
            huff_maps.append(huff_map)
            bitstrings.append(self._encode(inds, extra, huff_map))
        
        self._bitstrings = bitstrings
        self._huff_maps = huff_maps
        
    def decompress(self):
        groups = []
        for i in xrange(3):
            indices, extras = self._decode(self._bitstrings[i], self._huff_maps[i])
            groups.append(self._indicesToCoeffs(indices, extras))
            
        q_subbands = self._ungroup(groups, self._shapes, self._tvals)
        
        subbands = [self._dequantize(q_subbands[i], self._Q[i], self._Z[i]) for i in xrange(64)]
        
        img = self._recreate(subbands)
        return self._postProcess()
    
    def getRatio(self):
        bitlength = len(self._bitstrings[0]) + len(self._bitstrings[1]) + len(self._bitstrings[2])
        return float(self._pixels*8) / float(bitlength)

finger = imread('finger.pgm')
def compare(r, img=finger):
    wsq = WSQ()
    wsq.compress(finger, r)
    print 'For r of ' + str(r) + ', we get a compression ratio of ' + str(wsq.getRatio())
    new_finger = wsq.decompress()
    plt.subplot(211)
    plt.imshow(finger, cmap=plt.cm.Greys_r)
    plt.subplot(212)
    plt.imshow(new_finger.clip(finger.min(), finger.max()), cmap=plt.cm.Greys_r)
    plt.show()

def make_fig():
    wsq1 = WSQ()
    wsq1.compress(finger,0.5)
    new1 = wsq1.decompress()

    wsq2 = WSQ()
    wsq2.compress(finger,0.25)
    new2 = wsq2.decompress()

    wsq3 = WSQ()
    wsq3.compress(finger,0.12)
    new3 = wsq3.decompress()

    wsq4 = WSQ()
    wsq4.compress(finger,0.05)
    new4 = wsq4.decompress()

    wsq5 = WSQ()
    wsq5.compress(finger,0.025)
    new5 = wsq5.decompress()

    plt.subplot(231)
    plt.xticks(np.array([]))
    plt.yticks(np.array([]))
    plt.xlabel('Original Image')
    plt.imshow(finger, cmap=plt.cm.Greys_r)

    plt.subplot(232)
    plt.xticks(np.array([]))
    plt.yticks(np.array([]))
    plt.xlabel('Compression Ratio 20:1')
    plt.imshow(new1.clip(finger.min(), finger.max()), cmap=plt.cm.Greys_r)

    plt.subplot(233)
    plt.xticks(np.array([]))
    plt.yticks(np.array([]))
    plt.xlabel('Compression Ratio 45:1')
    plt.imshow(new2.clip(finger.min(), finger.max()), cmap=plt.cm.Greys_r)
    
    plt.subplot(234)
    plt.xticks(np.array([]))
    plt.yticks(np.array([]))
    plt.xlabel('Compression Ratio 100:1')
    plt.imshow(new3.clip(finger.min(), finger.max()), cmap=plt.cm.Greys_r)
    
    plt.subplot(235)
    plt.xticks(np.array([]))
    plt.yticks(np.array([]))
    plt.xlabel('Compression Ratio 223:1')
    plt.imshow(new4.clip(finger.min(), finger.max()), cmap=plt.cm.Greys_r)
    
    plt.subplot(236)
    plt.xticks(np.array([]))
    plt.yticks(np.array([]))
    plt.xlabel('Compression Ratio 450:1')
    plt.imshow(new5.clip(finger.min(), finger.max()), cmap=plt.cm.Greys_r)
    plt.savefig('SRC_fingerprints.jpg')