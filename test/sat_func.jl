using Test
using  TensorBranching: cal_gamma_0, cal_gamma_1, cal_gamma_2, get01s

n=30
m=100
times=1e3
t=0
max_g=zeros(1,3)
sum_g=zeros(1,3)

#核心算法是g0，g1是兜底，g2是防止极端情况

#下面进行times次随机模拟，计算三种算法的最差和平均gamma
#注意：由于"SAT_func.jl"中随机生成的01串之间需要互不相同，因此测试时不应让m接近2^n, 否则随机生成strs的函数效率会非常低
while t<times
    strs=get01s(n,m)
    strs=trunc.(Int64,strs)
    g1=cal_gamma_1(n,m)
    g2=cal_gamma_2(n,m,strs)
    g0=cal_gamma_0(n,m,strs)
    max_g=[max(max_g[1],g1) max(max_g[2],g2) max(max_g[3],g0)]
    sum_g+=[g1 g2 g0]
    t+=1
end
println(max_g)
println(sum_g/times)
    
# 下面的代码是为了证明，即使出现诸如strs={e_i}(e_i是只有第i个位置为1，其余位置为0的向量)，cal0依然优于cal1和cal2
N=20
testA=zeros(N,N)
for i=1:N
    testA[i,i]=1
end
cal_gamma_1(N,N)
cal_gamma_2(N,N,testA)
cal_gamma_0(N,N,testA)