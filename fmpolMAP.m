function [mu,sig2,fms,M,CGE] = fmpolMAP(C0,D0,alpha,T,kms,R,optvct)

ss=size(C0); p=ss(1); sss=size(ss); 
Tsteps=1; if sss(2)>2, Tsteps=ss(3); end; 
ss=size(kms); ks=ss(2); kmax=max(kms); km=kmax+1; km3=max(3,km); 
F=zeros(p*km3,p*km3); M=zeros(p,p,km3); fms=zeros(1,km3); nn=zeros(1,R);

H=eye(p*km3); jrepls=randsample(p,R,true,alpha); 
for Ts=1:Tsteps
    if Tsteps>1, C=C0(:,:,Ts);  D=D0(:,:,Ts); else C=C0; D=D0; end;
    for k=1:km3
        F((k-1)*p+1:k*p,(k-1)*p+1:k*p)=C+D;
        if k<km3, F((k-1)*p+1:k*p,k*p+1:(k+1)*p)=D; end;
    end;
    H=H*expm(F*T/Tsteps);
    for repl=1:R, 
        [n,jrepls(repl)]= MAPsim(p,C,D,jrepls(repl),T/Tsteps);
        nn(1,repl)=nn(1,repl)+n; if mod(repl,10^4)==0, [Ts repl], end;
    end;
end;


for k=1:km3
    Mk=H(1:p,(k-1)*p+1:k*p); M(:,:,k)=Mk;
    fms(k)=alpha*Mk*ones(p,1)*factorial(k-1);
end; 

mu=fms(2); sig2=fms(3)+fms(2)-fms(2)^2; 
fms=fms(1:km); muref=mu; sig2ref=sig2; figno=1;
xmax=ceil(1.2*(muref+2.58*sig2ref^(1/2))); poisson=mu>sig2;

if nargin>6
    ss=size(optvct); sopt=ss(2);
    if optvct(1)>0, figno=optvct(1); end;
    if sopt>1 & optvct(2)>0,  xmax=optvct(2); end;
    if sopt>2 & optvct(3)>0,  ymax=optvct(3); end;
    if sopt>3 & optvct(4)==1, poisson=1; muref=mu; end;
    if sopt>4 & optvct(5)==1, poisson=1; muref=sig2; end;
    if sopt>5 & optvct(6)>0,  poisson=1; muref=optvct(6); end;
    if sopt>6 & optvct(7)+optvct(8)>0
        poisson=0; muref=mu; sig2ref=sig2;
        if optvct(7)>0, muref=optvct(7); end;
        if sopt>7 & optvct(8)>0, sig2ref=optvct(8); end;
    end;
end;

xs=0:1:xmax; rrf=0*xs; frefx=0*xs; fsimx=0*xs;
CGE=zeros(ks,xmax+1);
for x=0:xmax, for repl=1:R
        if nn(repl)==x, fsimx(x+1)=fsimx(x+1)+1/R; end;
end; end;


if poisson, [A,c]=orthopolPC(muref,kmax,fms);
    frefx=exp(-muref)*muref.^xs./factorial(xs); 
else [r,rho,A,c]=orthopolNB(muref,sig2ref,kmax,fms); 
    frefx=(1-rho)^r*gamma(xs+r)/gamma(r)./factorial(xs).*rho.^xs;
end;

figure(figno); s=0; st={'-.k','--k','-k'};
if R>0, plot(xs,fsimx,'*r','LineWidth',4); hold on; end;
plot(xs,frefx,'--b','LineWidth',1); hold on;

for ik=1:ks 
    k=kms(ik); 
    [pdf,cdf]=expansion(k,A,c(1:k),xs,frefx); CGE(ik,:)=pdf;
    plot(xs,pdf,st{ik},'LineWidth',2); hold on;
end;
ymax=1.2*max(max(frefx),max(fsimx)); 
if nargin>6, if sopt>2 & optvct(3)>0, ymax=optvct(3); end; end;
ymin=-0.15*ymax; 
plot(xs,0*xs,'g','LineWidth',2);
axis([0 xmax ymin ymax]);
xleg=0.68*xmax:0.5:0.75*xmax; ysimleg=0*xleg+0.94*ymax; 

xsimleg=[0.69*xmax,0.73*xmax]; ysimleg=0*xsimleg+0.94*ymax;
plot(xsimleg,ysimleg,'*r','LineWidth',4);
if R>0, text(0.77*xmax,0.94*ymax,'Simulated','FontSize',18); end;
plot(xleg,0*xleg+0.88*ymax,'--b','LineWidth',1);
text(0.77*xmax,0.88*ymax,'Reference','FontSize',18);
for ik=1:ks
    plot(xleg,0*xleg+(0.88-0.06*ik)*ymax,st{ik},'LineWidth',2);
    text(0.77*xmax,(0.88-0.06*ik)*ymax,'N = ','FontSize',18);
    text(0.85*xmax,(0.88-0.06*ik)*ymax,num2str(kms(ik)),'FontSize',18);
end;
if xmax<25, xticks(0:5:35); else xticks('auto'); end;
%xticklabels({'x=0','x=5','x=10'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,c] = orthopolPC(lambda,kmax,fms)
km=kmax+1; c=0*(1:km); c(1)=1;
A=zeros(km,km); A(1,1)=1;
for n=2:km
    for k=1:n
        A(n,k)=(-1)^(n-k)*lambda^((n-1)/2-k+1);
        A(n,k)=A(n,k)/factorial(k-1)/factorial(n-k);
        A(n,k)=A(n,k)*sqrt(factorial(n-1));
    end;
    s=0; for k=1:n, s=s+A(n,k)*fms(k); end; c(n)=s;
end;

end %ortopolPC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r,rho,A,c] = orthopolNB(mu,sig2,kmax,fms)
km=kmax+1; c=0*(1:km); c(1)=1;
A=zeros(km,km); A(1,1)=1;
rho=1-mu/sig2; r=T*(1-rho)/rho;
A=zeros(km,km);  A(1,1)=1; 
raf=0*(1:km); raf(1)=1; raf(2)=r; % rising factorial
for k=3:km, raf(k)=raf(k-1)*(r+k-2); end;
ndf=zeros(km,km); % falling factorial
ndf(1,1)=1;
for n=2:km 
    ndf(n,1)=1; ndf(n,2)=n-1;
    for k=3:n, ndf(n,k)=ndf(n,k-1)*(n-k+1); end;
end;
for n=1:km
    for k=1:n
        A(n,k)=(rho^(n-1)*raf(n)/factorial(n-1))^(1/2);
        A(n,k)=A(n,k)*(-1)^(k-1)*ndf(n,k)*(1-rho)^(k-1);
        A(n,k)=A(n,k)/raf(k)/rho^(k-1)/factorial(k-1);
    end;
    s=0; for k=1:n, s=s+A(n,k)*fms(k); end; c(n)=s;
end;
end %ortopol NB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pdf,cdf] = expansion(km,A,c,t,f0t)
%tk(1,:)=1, tk(k,:)=t(t-1)...(t-k+2) for k>1
nn=size(t); N=nn(2); km=km;
tk=zeros(km,N); tk(1,:)=t*0+1;
for k=2:km, tk(k,:)=tk(k-1,:).*(t-k+2); end; 
incr=c(1)*tk(1,:); 
for n=2:km
    incrn=0*t;
    for k=1:n, incrn=incrn+A(n,k)*tk(k,:); end; 
    incr=incr+c(n)*incrn;
end;
pdf=f0t.*incr; cdf=cumsum(incr);
end %expansion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n,j]=MAPsim(p,C,D,j0,T)
pint=1.2*max(abs(diag(C))); 
probs=zeros(p,2*p+1);
for i=1:p
    for k=1:p
        if k~=i, probs(i,k)=probs(i,k)+C(i,k); end;
        probs(i,p+k)=probs(i,p+k)+D(i,k);
    end; 
    probs(i,2*p+1)=pint-sum(probs(i,:)); 
end;
probs=probs/pint; 
n=0; t=0; j=j0;
while t<T
    t=t+exprnd(1)/pint;
    if t<T
        k=randsample(2*p+1,1,true,probs(j,:));
        if k<=p, j=k; elseif k<=2*p, j=k-p; n=n+1; end;
    end;
end;
%if (k-p)*(2*p+1-k)>0, n=n-1; end;
end %MAPsim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %fmpolMAP

