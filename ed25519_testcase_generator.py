import random
q = pow(2,255)-19
q_inv = 21330121701610878104342023554231983025602365596302209165163239159352418617883 # q*q_inv % 2^255 = 2^255-1 = -1 mod 2^255
R = pow(2,255)%q
count_add_sub = 0
count_mul = 0

class number:
    def __init__(self, value: int):
        self.value = value # 255 bit
        
    def __add__(self, other: 'number') -> 'number': 
        r = self.value + other.value
        if(r>=q):
            r -= q
        assert r == ((self.value + other.value) % q)
        global count_add_sub
        count_add_sub += 1
        return number(r)

    def __sub__(self, other: 'number') -> 'number':
        if(self.value >= other.value):
            r = self.value - other.value
        else:
            r = q - other.value
            r += self.value
        assert r == ((self.value - other.value) % q)
        global count_add_sub
        count_add_sub += 1
        return number(r)
    
    def MM(self,value1: int, value2: int) -> int: # Montgomery multiplication: (value1 * value2)>>255 mod q
        r = value1 * value2
        tmp = (((r%pow(2,255))*q_inv)%pow(2,255))*q
        r = (r + tmp)>>255
        if(r>=q):
            r -= q
        global count_mul
        count_mul += 1
        return r
    
    def __mul__(self, other: 'number') -> 'number': # mod mul: value1 * value2 mod q
        r = self.MM(self.value,R*R%q)
        r = self.MM(r,other.value)
        assert r == ((self.value * other.value) % q)
        return number(r)
    
    def __truediv__(self, other: 'number') -> 'number': # mod div: value1 / value2 mod q
        #calculate value2 ^ (p-2) mod p
        q_minus_2_in_bin = "111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111101011"
        r = number(1)
        for i in range(255):
            r = r*r
            if(q_minus_2_in_bin[i]=="1"):
                r = r*other
        #calculate value1/value2 mod p
        r = self*r
        return r
    
    def __eq__(self, other: 'number') -> bool: 	# used for debug
        return self.value==other.value
	

d = number(0x52036cee2b6ffe738cc740797779e89800700a4d4141d8ab75eb4dca135978a3)
class point:
	def __init__(self, number_X: number, number_Y: number, number_Z: number = number(1)):
		self.X = number_X
		self.Y = number_Y
		self.Z = number_Z

	def __add__(self, other: 'point') -> 'point':
		Z1Z2 = self.Z*other.Z
		X1X2Y1Y2 = self.X*other.X*self.Y*other.Y
		X3 = Z1Z2*(self.X*other.Y+other.X*self.Y)*(Z1Z2*Z1Z2-d*X1X2Y1Y2)
		Y3 = Z1Z2*(self.Y*other.Y+self.X*other.X)*(Z1Z2*Z1Z2+d*X1X2Y1Y2)
		Z3 = (Z1Z2*Z1Z2-d*X1X2Y1Y2)*(Z1Z2*Z1Z2+d*X1X2Y1Y2)
		return point(X3,Y3,Z3)

	def __mul__(self, M: int) -> 'point':
		r = point(number(0), number(1))  # the zero point
		M_in_bin = "{:0255b}".format(M)
		for i in range(255):
			r = r + r
			if(M_in_bin[i]=="1"):
				r = r + self
		return r

	def reduce(self) -> 'point':
		x = self.X/self.Z
		y = self.Y/self.Z
		if(x.value%2==1): x.value = q-x.value
		if(y.value%2==1): y.value = q-y.value
		return point(x, y)
		
	def is_on_curve(self): # require self.Z=1
		return self.Y*self.Y-self.X*self.X == number(1) + d * self.X*self.X*self.Y*self.Y
	
	def __str__(self): # used for debug
		if(self.Z.value!=1):
			text = "The z-coordinate != 1"
		elif(self.is_on_curve()):
			text = "x: {:064x}\n".format(self.X.value) + "y: {:064x}\n".format(self.Y.value)
		else:
			text = "Invalid point"
		return text


# Modular square root given by GPT
from sympy.ntheory import legendre_symbol
def tonelli_shanks(k, p):
    """
    计算 k 的平方根 (mod p)，假设 p 是质数。
    通过 Tonelli-Shanks 算法来实现。
    """
    # 如果 k ≡ 0 (mod p)，那么 y 也应该是 0
    if k == 0:
        return 0
    
    # 如果 p ≡ 3 (mod 4)，可以使用快速方法来直接求解
    if p % 4 == 3:
        return pow(k, (p + 1) // 4, p)
    
    # Tonelli-Shanks 算法
    s, q = 0, p - 1
    while q % 2 == 0:
        s += 1
        q //= 2
    
    # 找到一个非二次剩余的 z (mod p)
    z = 2
    while pow(z, (p - 1) // 2, p) == 1:
        z += 1
    
    # 初始化变量
    m = s
    c = pow(z, q, p)
    t = pow(k, q, p)
    r = pow(k, (q + 1) // 2, p)
    
    while t != 0 and t != 1:
        # 寻找 t 的最小非零平方根
        t2i = t
        i = 0
        for i in range(1, m):
            t2i = pow(t2i, 2, p)
            if t2i == 1:
                break
        
        # 更新 m、c、t、r
        b = pow(c, 2**(m - i - 1), p)
        m = i
        c = pow(b, 2, p)
        t = (t * b * b) % p
        r = (r * b) % p
    
    if t == 0:
        return 0
    else:
        return r


if __name__ == "__main__":
    n_case = 1
    
    #### Generate testcase ####
    M_list = []
    x_list = []
    y_list = []
    G_list = []
    i = 0
    # Get valid (x, y) with random x
    while i < n_case:
        x = random.randint(0, q-1)
        x_2 = number(x)*number(x)   # x^2
        y_2 = (number(1) - (number(q-1)*x_2)) / (number(1) - (d*x_2)) # y^2
        if legendre_symbol(y_2.value, q)==1:    # y^2 have integer square root
            i += 1
            M_list.append( random.randint(0, (1 << 255) - 1) )
            y = tonelli_shanks(y_2.value, q)
            x_list.append(x)
            y_list.append(y)

    # # Get valid (x, y) with random y
    # while i < n_case:
    #     y = random.randint(0, q-1)
    #     y_2 = number(y)*number(y)   # y^2
    #     x_2 = (number(1) - (y_2)) / (number(q-1) - (d*y_2)) # x^2
    #     if legendre_symbol(x_2.value, q)==1:    # x^2 have integer square root
    #         i += 1
    #         M_list.append( random.randint(0, (1 << 255) - 1) )
    #         x = tonelli_shanks(x_2.value, q)
    #         x_list.append(x)
    #         y_list.append(y)
            
    # Compute answer
    for i in range(n_case):
        point_P = point(number(x_list[i]), number(y_list[i]))
        point_G = (point_P * M_list[i]).reduce()
        G_list.append( [point_G.X.value, point_G.Y.value] )
        
    #### Write file ####
    with open("tb_dat_h.sv", "w", encoding="utf-8") as file:
        file.writelines("package dat_h;\n")
        file.writelines("integer pat_num = 3;\n")
        file.writelines("reg [" +str(768*n_case-1)+ ":0] input_data  = " +str(768*n_case)+ "'h")
        for i in range(n_case):
            file.writelines(f"{M_list[i]:064x}{x_list[i]:064x}{y_list[i]:064x}")
        file.writelines(";\nreg [" +str(512*n_case-1)+ ":0] golden_data = " +str(512*n_case)+ "'h")
        for i in range(n_case):
            file.writelines(f"{G_list[i][0]:064x}{G_list[i][1]:064x}")
        file.writelines(";\nendpackage")
        
        
    # #### Add/Sub testcase ####
    # answer_add_list = []
    # answer_sub_list = []
    # for i in range(n_case):
    #     answer_add_list.append((number(x_list[i])+number(y_list[i])).value)
    #     answer_sub_list.append((number(x_list[i])-number(y_list[i])).value)
    # with open("add_sub_pattern.sv", "w", encoding="utf-8") as file:
    #     file.writelines("package dat_0;\n")
    #     file.writelines("integer pat_num = 0;\n")
    #     file.writelines("reg [" +str(256*n_case-1)+ ":0] data1 = " +str(256*n_case)+ "'h")
    #     for i in range(n_case):
    #         file.writelines(f"{x_list[i]:064x}")
    #     file.writelines(";\nreg [" +str(256*n_case-1)+ ":0] data2 = " +str(256*n_case)+ "'h")
    #     for i in range(n_case):
    #         file.writelines(f"{y_list[i]:064x}")
    #     file.writelines(";\nreg [" +str(256*n_case-1)+ ":0] add_answer = " +str(256*n_case)+ "'h")
    #     for i in range(n_case):
    #         file.writelines(f"{answer_add_list[i]:064x}")
    #     file.writelines(";\nreg [" +str(256*n_case-1)+ ":0] sub_answer = " +str(256*n_case)+ "'h")
    #     for i in range(n_case):
    #         file.writelines(f"{answer_sub_list[i]:064x}")   
    #     file.writelines(";\nendpackage")
    
    # #### MM testcase ####
    # answer_mm_list = []
    # for i in range(n_case):
    #     answer_mm_list.append((number(x_list[i])*number(y_list[i])).value)
    # with open("mm_pattern.sv", "w", encoding="utf-8") as file:
    #     file.writelines("package dat_0;\n")
    #     file.writelines("integer pat_num = 0;\n")
    #     file.writelines("reg [" +str(256*n_case-1)+ ":0] data1 = " +str(256*n_case)+ "'h")
    #     for i in range(n_case):
    #         file.writelines(f"{x_list[i]:064x}")
    #     file.writelines(";\nreg [" +str(256*n_case-1)+ ":0] data2 = " +str(256*n_case)+ "'h")
    #     for i in range(n_case):
    #         file.writelines(f"{y_list[i]:064x}")
    #     file.writelines(";\nreg [" +str(256*n_case-1)+ ":0] mm_answer = " +str(256*n_case)+ "'h")
    #     for i in range(n_case):
    #         file.writelines(f"{answer_mm_list[i]:064x}")
    #     file.writelines(";\nendpackage")

        
        

