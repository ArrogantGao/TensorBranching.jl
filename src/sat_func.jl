
# SATBranch
# To optimize SAT by branching
# It includes some algorithm to optimize SAT by braching. Codes to relize these algorithm are in SAT_func.jl.

# To solve a SAT problem of size $N$ we need $\gamma^N$, what we do is to design braching algorithm whose $\gamma$ is as small as possible. 

# Now there is 3 algorithm in SAT_func.jl, in which cal_gamma_1 and cal_gamma_2 are trivial, as a reference. Important algorithm is cal_gamma_0, getting small $\gamma$ by Merge similar 01 strings.

# Example gives some test of these algorithms.
#The function soleq is used to solve equations 1=x^{-d1}+..+x^{-dk} using binary methods
function soleq(poly)
    k=size(poly,2)
    eps=1e-9
    low=1
    high=2
    while high-low>eps
        mid=(high+low)/2
        poly_res=0
        for i=1:k
            if poly[i]<=0 
                continue
            end
            poly_res+=mid^(-poly[i])
        end
        if poly_res>1
            low=mid
        else
            high=mid
        end
    end
    return (low+high)/2
end

#Randomly generate m～[30,120] and m different 01 strings
function get01s(n,m)
    if m>2^n
        println("error")
        return 
    end
    strs=zeros(m,n)
    i=1
    while i<=m
        s=rand((0,1),n);
        flag=1
        for j=1:i-1
            if s==strs[j,:]
                flag=0
            end
        end
        if flag==1
            strs[i,:]=s
            i+=1
        end
    end
    return strs
end

#合并01串的思路
function dist(str1,str2)
    s=0
    size_str=0
    for i in str1
        size_str+=1
    end
    for i=1:size_str
        if str1[i]<-0.5||str2[i]<-0.5
            s+=1
            continue
        end
        s+=abs(str1[i]-str2[i])
    end
    return s
end

function cal_gamma_0(n,m,strs)
    poly=ones(1,m)*n
    gamma=soleq(poly)
    d=zeros(m,m);
    st=trunc.(Int64,ones(m))
    #接下里对海明距离最小的两个01串进行合并尝试，若gamma减小了，合并并删除其中一个

    #先求两两之间的海明距离矩阵d
    for i=1:m
        for j=1:m
            d[i,j]=dist(strs[i,:],strs[j,:])
        end
    end
    for i=1:m
        d[i,i]=n+10
    end

    #接下来开始迭代，每次找出距离最小的两点，记录下所有的gamma
    while true
        mx=my=0
        mv=n+1
        for i=1:m
            for j=1:m
                if st[i]==1&&st[j]==1
                    if d[i,j]<mv
                        mv=d[i,j]
                        mx,my=i,j
                    end
                end
            end
        end
        if mv>=n
            break
        end
        #接下来开始处理poly
        poly[my]=0
        poly[mx]=mv
        gamma=min(soleq(poly),gamma)
        st[my]=0
        for i=1:n
            if strs[mx,i]!=strs[my,i]
                strs[mx,i]=-1
            end
        end
        for i=1:m
            if st[i]==1&&i!=mx
                d[i,mx]=d[mx,i]=dist(strs[mx,:],strs[i,:])
            end
        end
    end
    return gamma
end

#m*clause(1,n)算法
function cal_gamma_1(n,m)
    return m^(1/n)
end

#clause(1,n)+clause(m-1,1)
function cal_gamma_2(n,m,strs)
    strs=2*strs.-1
    mv=findmin(abs.(sum(strs,dims=1)))[1]
    mv=trunc(Int64,(m-mv)/2)

    poly=[ones(1,mv)*n 1] 
    return soleq(poly)
end







