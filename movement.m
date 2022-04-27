N_Sites=20;%site number
m=1;%mass
Dim_Extended=2^N_Sites; 
TimeSavingMode=false; % Use existing variables. 
insertion_site=10;

%Hamiltonian
H0=1/2*SpinChainH('X',N_Sites,'pbc')-m*SpinChainH('Z',N_Sites,'pbc');
if (~exist('IsPhysical','var'))||(~TimeSavingMode)
IsPhysical=true(Dim_Extended,1);
for i_Site=1:N_Sites
    Projector=~diag(((speye(Dim_Extended)+OnsitePauliMString('Z',i_Site,N_Sites))...
        .*(speye(Dim_Extended)+OnsitePauliMString('Z',mod(i_Site,N_Sites)+1,N_Sites))));
    IsPhysical=IsPhysical&Projector;
end
end
Dim_Physical=full(sum(IsPhysical));

H_eff=H0(IsPhysical,IsPhysical);
[groundstate, genergy] = eigs(H_eff,1,'smallestreal');

n1=-1/2*(OnsitePauliMString(3,insertion_site,N_Sites)+OnsitePauliMString(3,insertion_site+1,N_Sites));
%p1=sum(conj(groundstate).*(n1(IsPhysical,IsPhysical)*groundstate));

%s2=ones(N_Sites,1);
%s2(4:2:N_Sites)=0;
%s2_index=BasisIndex(s2,IsPhysical);
%p2=abs(groundstate(s2_index)).^2;
%p2/p1

if m<0
    psi=n1(IsPhysical,IsPhysical)*groundstate;
    v0=zeros(2^N_Sites,1);
    v0(IsPhysical)=1/sqrt(psi'*psi)*psi;

else
    psi=(eye(Dim_Physical)-n1(IsPhysical,IsPhysical))*groundstate;
    v0=zeros(2^N_Sites,1);
    v0(IsPhysical)=1/sqrt(psi'*psi)*psi;
end

if m<0
  L1=true(2^N_Sites,1);
    for i=1:N_Sites
        if i==insertion_site
            F=5/2+diag(1/2*OnsitePauliMString('ZZ',i,N_Sites)+3/2*OnsitePauliMString('Z',i,N_Sites)+3/2*OnsitePauliMString('Z',mod(i,N_Sites)+1,N_Sites));
        else
            F=1+diag(OnsitePauliMString('ZZ',i,N_Sites)+OnsitePauliMString('Z',i,N_Sites)+OnsitePauliMString('Z',mod(i,N_Sites)+1,N_Sites));
        end
        L1=L1&(F==0);
    end
else
  L1=true(2^N_Sites,1);
    for i=1:N_Sites
        if i==insertion_site
            F=1/2+diag(1/2*OnsitePauliMString('ZZ',i,N_Sites)-1/2*OnsitePauliMString('Z',i,N_Sites)-1/2*OnsitePauliMString('Z',mod(i,N_Sites)+1,N_Sites));
        else
            F=1+diag(OnsitePauliMString('ZZ',i,N_Sites)+OnsitePauliMString('Z',i,N_Sites)+OnsitePauliMString('Z',mod(i,N_Sites)+1,N_Sites));
        end
        L1=L1&(F==0);
    end  
end

    H1=H0(L1,L1);
    N_eigenval=100;
    [ev, eigenenergy] = eigs(H1,N_eigenval,'smallestreal');
    v1=zeros(2^N_Sites,N_eigenval);
    v1(L1,:)=ev;

    t_start = 0;
    t_end = 50;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ;
    step = .1;
    y = t_start:step:(t_end-t_start-step);

    psi=v1*(exp(-1i*diag(eigenenergy)*y).*(v1'*v0));

    charge=zeros(N_Sites,length(y));

    for j=1:N_Sites

        charge(j,:)=1/2*sum(conj(psi).*((-OnsitePauliMString(3,mod(j,N_Sites)+1,N_Sites)-OnsitePauliMString(3,j,N_Sites))*psi),1);

    end

    test_charge=zeros(N_Sites,1);
    if m<0
        test_charge(insertion_site)=-1;
    else
        test_charge(insertion_site)=1;
    end

    phys_charge=charge+test_charge;
    h=surf(1:N_Sites,y,phys_charge.');
    set(h,'edgecolor','none');
    view([0,0,1]);
    set(gcf,'Renderer','Painters');
    xlabel('$n$','fontsize',14,'Interpreter','latex');
    ylabel('$t$','fontsize',14,'Interpreter','latex');
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'fontsize',14);
    colorbar;
    caxis([0,1]);
    TheColorBar=colorbar;
    TheColorBar.Label.String='charge';
    TheColorBar.Label.Interpreter='latex'; 
    TheColorBar.TickLabelInterpreter='latex';
    xlim([1,N_Sites]);
    