N=20;%site number
m=-0.1;
insertion_site=10;

H0=1/2*SpinChainH('X',N,'pbc')-m*SpinChainH('Z',N,'pbc');
L0=true(2^N,1);
    for i=1:N
        F=1+diag(OnsitePauliMString('ZZ',i,N)+OnsitePauliMString('Z',i,N)+OnsitePauliMString('Z',mod(i,N)+1,N));
        L0=L0&(F==0);
    end

H_eff=H0(L0,L0);
[groundstate, genergy] = eigs(H_eff,1,'smallestreal');
[eigenstate, eigenenergy] = eigs(H_eff,200,'smallestreal');
sigmax=OnsitePauliMString('X',insertion_site,N);
identity=OnsitePauliMString('I',insertion_site,N);
pip=1/4*(OnsitePauliMString('I',insertion_site,N)-OnsitePauliMString('Z',(insertion_site-1),N)-OnsitePauliMString('Z',(insertion_site+1),N)+OnsitePauliMString('ZIZ',(insertion_site-1),N));
initialstate=(sigmax(L0,L0)+identity(L0,L0)-pip(L0,L0))*groundstate;

t_start = 0;
t_end = 25;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ;
step = .05;
y = t_start:step:(t_end-t_start-step);

psi=eigenstate*(exp(-1i*diag(eigenenergy)*y).*(eigenstate'*initialstate));

spin=zeros(N,length(y));
nb=zeros(N,1);
sping=zeros(N,1);
spini=zeros(N,1);

for j=1:N
    sigmaz=OnsitePauliMString(3,j,N);
    noperator=OnsitePauliMString(3,j,N)+OnsitePauliMString(3,mod(j,N)+1,N);
    spin(j,:)=1/2*sum(conj(psi).*(sigmaz(L0,L0)*psi),1);
    sping(j)=1/2*sum(conj(groundstate).*(sigmaz(L0,L0)*groundstate),1);
    spini(j)=1/2*sum(conj(initialstate).*(sigmaz(L0,L0)*initialstate),1);
    nb(j)=-1/2*sum(conj(groundstate).*(noperator(L0,L0)*groundstate),1);

end

 nbar=nb*ones(1,length(y));
 c=zeros(N,N/2+2,length(y));
 
  
 for p=1:(N/2+2)
     j=p-1;
     for k=1:N
         nk=-1/2*(OnsitePauliMString(3,k,N)+OnsitePauliMString(3,mod(k,N)+1,N));
         nkj=-1/2*(OnsitePauliMString(3,mod((k+j)-1,N)+1,N)+OnsitePauliMString(3,mod((k+j),N)+1,N));
         c(k,p,:)=sum(conj(psi).*(nk(L0,L0)*nkj(L0,L0)*psi),1)-sum(conj(psi).*(nk(L0,L0)*psi),1).*nbar(mod((j+k-1),N)+1,:)-sum(conj(psi).*(nkj(L0,L0)*psi),1).*nbar(k,:)+nbar(mod((j+k-1),N)+1,:).*nbar(k,:);
     end
 
 end
 
% c1=zeros(N,N/2);
% c2=zeros(N,N/2);
% c3=zeros(N,N/2);
% c4=zeros(N,N/2);
% c5=zeros(N,N/2);
% for j=1:(N/2)
%     for k=1:N
%         nk=-1/2*(OnsitePauliMString(3,k,N)+OnsitePauliMString(3,mod(k,N)+1,N));
%         nkj=-1/2*(OnsitePauliMString(3,mod((k+j)-1,N)+1,N)+OnsitePauliMString(3,mod((k+j),N)+1,N));
%         c1(k,j)=sum(conj(initialstate).*(nk(L0,L0)*nkj(L0,L0)*initialstate),1);
%         c2(k,j)=sum(conj(initialstate).*(nk(L0,L0)*initialstate),1).*nb(mod((j+k-1),N)+1);
%         c3(k,j)=sum(conj(initialstate).*(nkj(L0,L0)*initialstate),1).*nb(k);
%         c4(k,j)=nb(mod((j+k-1),N)+1).*nb(k);
%         c5(k,j)=sum(conj(initialstate).*(nk(L0,L0)*initialstate),1)*sum(conj(initialstate).*(nkj(L0,L0)*initialstate),1);
%     end
% 
% end



 f=reshape(sum((c),1),(N/2+2),length(y));
  h=surf(0:(N/2+1),y,f');
      set(h,'edgecolor','none');
      view([0,0,1]);
      set(gcf,'Renderer','Painters');
      xlabel('$r$','fontsize',14,'Interpreter','latex');
      ylabel('$t$','fontsize',14,'Interpreter','latex');
      set(gca,'TickLabelInterpreter', 'latex');
      set(gca,'fontsize',14);
      colorbar;
      TheColorBar=colorbar;
      TheColorBar.Label.Interpreter='latex'; 
      TheColorBar.Ticks;
      xlim([0,N/2+1]);
      

