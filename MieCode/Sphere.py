def func():
	print("hi Kevin")

class Elements:
    elementCount = 0
    
    def __init__(self, chemSymb, wp, gamma, vf, A, RefracFile):
        self.chemSymb = chemSymb
        self.wp = wp #E15 Hz
        self.gamma = gamma #s
        self.vf = vf #m/s
        self.A = A
        self.RefracFile = RefracFile
        Elements.elementCount += 1