from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, QISKitError, QuantumJob
from qiskit import available_backends, execute, register, get_backend, compile
from qiskit.tools.visualization import circuit_drawer
from pprint import pprint
from math import pi
import matplotlib.pyplot as plt
import time 
import json
import numpy as np

##########Some useful classical function##########

# pre_chwck function:
# to check if input n is even or a non-trivial perfect power
def pre_check(n):
	if(n%2==0): return True
	l = int(np.ceil(np.log2(n)))
	# print(l)
	for power in np.arange(l,1,-1):
		for num in np.arange(3,int(np.sqrt(n))+1,2):
			ans = num ** power
			# print(num,power,ans)
			if(ans == n): return True
			if(ans > n): break;
	return False

# gcd function: return greatest common divisor(gcd) of two numbers
def gcd(n,m):
	# print("GCD:  ",n,m)
	if(m == 0): return n
	return gcd(m,int(n%m))

# cfe function: calculate continued fraction expansion
def cfe(a,b): #a/b
	a = int(a)
	b = int(b)
	ans = []
	while(True):
		# print(a,b)
		q = int(a/b)
		ans.append(q)
		c = a
		a = b
		b = c-b*q
		if(a == 1): break;
	return ans

# print(cfe(123,100))

def con(a,b):
	exp = cfe(a,b)
	# print("exp======" ,exp)
	p = []
	q = []
	p.append(exp[0])
	p.append(exp[0]*exp[1]+1)
	q.append(1)
	q.append(exp[1])
	for i in np.arange(2,len(exp),1):
		p.append(exp[i]*p[i-1]+p[i-2])
		q.append(exp[i]*q[i-1]+q[i-2])
	ans = [p,q]
	return ans

# print(con(3,4))

def get_approx(s0,r0):
	p_list = con(s0,r0)[0]
	q_list = con(s0,r0)[1]
	for i in range(len(p_list)):
		if(q_list[i]%2==0 and not (2*q_list[i]**2*abs(r0*p_list[i]-s0*q_list[i]) > q_list[i]*r0)):
			return[p_list[i],q_list[i]]
	return [q_list[-1],p_list[-1]]



## N is the number we want to factorize
## x is the randomly chosen number
## s0/r0 is what we get from continued fraction expansion
def factor(N,x,s0,r0):
	r = int(get_approx(s0,r0)[1])
	if(r%2==0):
		a = int(gcd(x**(r/2)+1,N))
		b = int(gcd(x**(r/2)-1,N))
		if(a != 1): return a
		if(b != 1): return b
	return 1	


##########The core of QISkit##########

# Access ibmqx
qx_config = {
    "APItoken": YOUR_TOKEN,
    "url":"https://quantumexperience.ng.bluemix.net/api"}
try:
    register(qx_config['APItoken'], qx_config['url'])
    print('\nYou have access to great power!')
except: 
    print('Something went wrong.\nDid you enter a correct token?')

# Some parameter
x_list = {2, 7, 8, 11, 13};
data = {};
shot = 1024
reg1 = 3 #9
reg2 = 4
backend = "local_qasm_simulator" #"ibmqx5"

# Q circuit
for x in x_list:
    q = QuantumRegister(reg1 + reg2)
    c = ClassicalRegister(reg1 + reg2)
    qc = QuantumCircuit(q, c)
    
    # init
    qc.x(q[reg1 + reg2 - 1])
    
    # superposition
    for i in range(0, reg1):
        qc.h(q[i])
    
    # modular register 2
    for i in range(0, reg1):
        for j in range(0, 2 ** i):
            control = q[i]
            if x == 2:
                qc.cswap(control, q[reg1 + 0], q[reg1 + 1])
                qc.cswap(control, q[reg1 + 1], q[reg1 + 2])
                qc.cswap(control, q[reg1 + 2], q[reg1 + 3])
            if x == 7:
                qc.cswap(control, q[reg1 + 2], q[reg1 + 3])
                qc.cswap(control, q[reg1 + 1], q[reg1 + 2])
                qc.cswap(control, q[reg1 + 0], q[reg1 + 1])
                qc.cx(control, q[reg1 + 0])
                qc.cx(control, q[reg1 + 1])
                qc.cx(control, q[reg1 + 2])
                qc.cx(control, q[reg1 + 3])
            if x == 8:
                qc.cswap(control, q[reg1 + 2], q[reg1 + 3])
                qc.cswap(control, q[reg1 + 1], q[reg1 + 2])
                qc.cswap(control, q[reg1 + 2], q[reg1 + 3])
            if x == 11:
                qc.cswap(control, q[reg1 + 1], q[reg1 + 3])
                qc.cswap(control, q[reg1 + 0], q[reg1 + 2])
                qc.cx(control, q[reg1 + 0])
                qc.cx(control, q[reg1 + 1])
                qc.cx(control, q[reg1 + 2])
                qc.cx(control, q[reg1 + 3])
            if x == 13:
                qc.cswap(control, q[reg1 + 0], q[reg1 + 1])
                qc.cswap(control, q[reg1 + 1], q[reg1 + 2])
                qc.cswap(control, q[reg1 + 2], q[reg1 + 3])
                qc.cx(control, q[reg1 + 0])
                qc.cx(control, q[reg1 + 1])
                qc.cx(control, q[reg1 + 2])
                qc.cx(control, q[reg1 + 3])
    # measurement register 1
    for i in range(0, reg1):
        qc.measure(q[i], c[i]);
    
    # IQFT
    # swap
    for i in range(0, int(reg1 / 2)):
        qc.swap(q[i] , q[reg1 - 1 - i])
    
    for i in range(0, reg1):
        for j in range(0, i):
            qc.cu1(-2*pi / (2**(i-j+1)), q[i], q[j])
        qc.h(q[i])
    qc.barrier(q)

    # measurement register 2
    for i in range(reg1, reg1 + reg2):
        qc.measure(q[i], c[i]);
    
# Draw Q circuit
    im = circuit_drawer(qc)
    im.save("circuit-" + str(x) + ".jpeg")

# Excute
    job_exp = execute(qc, backend, shots=shot)
    
# Wait job done
    lapse = 0
    interval = 5
    while not job_exp.done:
        print('Status @ {} seconds'.format(interval * lapse))
        print(job_exp.status)
        print(x)
        time.sleep(interval)
        lapse += 1
    print(job_exp.status)
    
    print(job_exp.result().get_counts(qc))
    data[x] = job_exp.result().get_counts(qc);
    
f = open("./result.json","w")
json.dump(data, f);

##########Classical analysis##########
f = open("./result3.json",'r')
load = json.load(f)
data = {'2':{},'7':{},'8':{},'11':{},'13':{}}

# Parse results
kk = load.keys()
for ik in kk:
    item = load[ik]
    kkk = item.keys()
    for key in kkk:
        r1=0
        r2=0
        for i in range(4, 7):
            r1 *= 2;
            r1 += (int(key[i]))
        for i in range(0, 4,):
            r2 *= 2;
            r2 += (int(key[i]))
        if (r2 not in data[ik].keys()):
            data[ik][r2] = [[],[]];
        data[ik][r2][1].append(r1);
        data[ik][r2][0].append(item[key])

# Calculate result, Plot spectrum
#plt.switch_backend('agg');
cnt = 1
for x in data:
    fail = 0
    for u2 in data[x]:
        title = 'x: '+str(x)+', u2: '+str(u2)
        amp = data[x][u2][0]
        r1 = data[x][u2][1]
        plt.figure(cnt)
        plt.stem(r1, amp)
        plt.title(title)
        plt.ylim(0,max(amp)*1.05)
        # plt.show()
        plt.savefig(title)
        cnt+=1
        for qqq in range(len(r1)):
            if(data[x][u2][1][qqq]==0): 
                fail += data[x][u2][0][qqq]
                continue
            gcdd = gcd(data[x][u2][1][qqq],8)
            ans = factor(15,int(x),int(data[x][u2][1][qqq]/gcdd),int(8/gcdd))
            # print (ans)
            if(ans != 3 and ans != 5):
                fail += data[x][u2][0][qqq]
    print("x= ", x, ", fail rate = ", fail/1024)
