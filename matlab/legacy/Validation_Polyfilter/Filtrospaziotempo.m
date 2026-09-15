function  [varout,varout1,varout2,varout3,varout4,varout5,varout6]= Filtrospaziotempo(Thot,Xrad)


varin=Thot;
grdx=squeeze(Thot(1,:,1))';
grdy=(squeeze(Thot(:,1,1)));
grdz=squeeze(Thot(1,1,:));

dimx=length(grdx); dimy=length(grdy); dimz=length(grdz);

varout=varin;

[J,I,K]=size(varin);
rad=min([2,floor(I/2)-1,floor(J/2)-1,floor(K/2)-1]);  % neighbourhood radius
rad=Xrad;
if rad>0

    if rad<2
        % linear regression
        for i=1:dimx
            for j=1:dimy
                for k=1:dimz
                    % set the neighbourhood indices including edge effects
                    iin=(i-rad)*((i-rad)>0)+((i-rad)<=0);
                    ifin=(i+rad)*((i+rad)<=dimx)+dimx*((i+rad)>dimx);
                    jin=(j-rad)*((j-rad)>0)+((j-rad)<=0);
                    jfin=(j+rad)*((j+rad)<=dimy)+dimy*((j+rad)>dimy);
                    kin=(k-rad)*((k-rad)>0)+((k-rad)<=0);
                    kfin=(k+rad)*((k+rad)<=dimz)+dimz*((k+rad)>dimz);

                    % select values belonging to the neighbourhood
                    varloc=varin(jin:jfin,iin:ifin,kin:kfin);
                    varloc1=varloc(:);


                    % select the local system of coordinates for the neighbourhood
                    xloc=[iin:ifin]-i;
                    yloc=[jin:jfin]-j;
                    zloc=[kin:kfin]-k;
                    [XLOC,YLOC,ZLOC]=meshgrid(xloc,yloc,zloc);
                    x1=XLOC(:);
                    x2=YLOC(:);
                    x3=ZLOC(:);

                    % array of coefficients for the second order polynomial fit
                    MAT=[ones(size(x1)) x1 x2 x3 ]; % buffer variable containing terms of the polynomial

                    coeff=MAT\varloc1;

                    % regressed value
                    varout(j,i,k)=coeff(1);
                end
            end
        end
    elseif rad>=2
        % quadratic regression
        fprintf('Quadratic regression. "i" value (1 --> 128):\n')
        for i=1:dimx
            fprintf('%d, ', i)
            for j=1:dimy
                for k=1:dimz
                    % set the neighbourhood indices including edge effects
                    iin=(i-rad)*((i-rad)>0)+((i-rad)<=0);
                    ifin=(i+rad)*((i+rad)<=dimx)+dimx*((i+rad)>dimx);
                    jin=(j-rad)*((j-rad)>0)+((j-rad)<=0);
                    jfin=(j+rad)*((j+rad)<=dimy)+dimy*((j+rad)>dimy);
                    kin=(k-rad)*((k-rad)>0)+((k-rad)<=0);
                    kfin=(k+rad)*((k+rad)<=dimz)+dimz*((k+rad)>dimz);

                    % select values belonging to the neighbourhood
                    varloc=varin(jin:jfin,iin:ifin,kin:kfin);
                    varloc1=varloc(:);


                    % select the local system of coordinates for the neighbourhood
                    xloc=[iin:ifin]-i;
                    yloc=[jin:jfin]-j;
                    zloc=[kin:kfin]-k;
                    [XLOC,YLOC,ZLOC]=meshgrid(xloc,yloc,zloc);
                    x1=XLOC(:);
                    x2=YLOC(:);
                    x3=ZLOC(:);

                    % array of coefficients for the second order polynomial fit
                    MAT=[ones(size(x1))...
                        x1 x2 x3 ...
                        x1.*x1 x2.*x2 x3.*x3 ...
                        x1.*x2 x1.*x3 x2.*x3 ]; % buffer variable containing terms of the polynomial

                    coeff=MAT\varloc1;

                    % regressed value
                    varout(j,i,k)=coeff(1);
                    varout1(j,i,k)=coeff(2);
                    varout2(j,i,k)=coeff(3);
                    varout3(j,i,k)=coeff(4);
                    varout4(j,i,k)=coeff(5);
                    varout5(j,i,k)=coeff(6);
                    varout6(j,i,k)=coeff(7);
                end
         
            end
        end

        
    end
else
    disp('      interrogation grid too small: ')
    disp('      regression not possible')
    varout=varin;
end


end