from scipy.fftpack import fft2, ifft2
from scipy.fftpack import fftshift, ifftshift
import matplotlib.pyplot as plt
from scipy import ndimage as nd
import numpy as np

image = 'einstein.jpeg'

im = nd.imread(image)
im = im[:,:,1]

'''
Plot powerspectrum
'''
'''
B = fft2(im)
B_shift = fftshift(B)
B_mag = np.abs(B_shift)
amp = 1000
B_mag *= amp/B_mag.max()
B_mag[B_mag>1] = 1

plt.subplot(121)
plt.imshow(im, cmap='gray')
plt.subplot(122)
plt.imshow(B_mag, cmap='gray')
plt.show()
'''

'''
Remove coefficients outside of a circle
'''

def fourier_radius(r, im):
	centerx = im.shape[0]/2
	centery = im.shape[1]/2
	#r = 100

	B = fft2(im)
	B_shift = fftshift(B)

	for x in xrange(im.shape[0]):
	    for y in xrange(im.shape[1]):
		if (x - centerx)**2 + (y - centery)**2 < r**2:
		    B_shift[x,y] = 0

	B_unshift = ifftshift(B_shift)
	new_im = ifft2(B_unshift)
	return new_im

plt.subplot(221)
plt.xticks(np.array([]))
plt.yticks(np.array([]))
plt.imshow(im, cmap='gray')
plt.subplot(222)
plt.xticks(np.array([]))
plt.yticks(np.array([]))
plt.imshow(np.abs(im - fourier_radius(100, im)), cmap='gray')
plt.subplot(223)
plt.xticks(np.array([]))
plt.yticks(np.array([]))
plt.imshow(np.abs(im - fourier_radius(50, im)), cmap='gray')
plt.subplot(224)
plt.xticks(np.array([]))
plt.yticks(np.array([]))
plt.imshow(np.abs(im - fourier_radius(25, im)), cmap='gray')
plt.tight_layout()
plt.savefig('SRC_einstein.jpg')






