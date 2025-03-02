%% All parameters are identical to those in Paper "Stacked Intelligent Metasurfaces for Multiuser Beamformingin the Wave Domain", whose code is available at https://github.com/JianchengAn/SIM-2-MUMISO.
close all
clear
clc
tic
%% System setting
Enlarge=0;
Bandwidth=2*10^(6); 
sigma=10^((-104)/10)*10^(2*Enlarge); 
Ferqun=28*10^9;
lambda=3*10^8/Ferqun;
MonteCarlo = 100;
Sum_Maxmin_rate=0;
%% User
User_Num=4; 
Coordinate_X_User=[17.2896339269092 29.9821853994423 43.5767305029477 57.4363251029044];%用户的横坐标
Coordinate_Y_User=[0 0 0 0 0 0];
Coordinate_Z_User=[0 0 0 0 0 0 ];
%% AP
Antenna_Num=User_Num;
d_antenna_spacing = 0.5*lambda; %
% Coordinate_X_AP=[0];
% Coordinate_Y_AP=[-30];
% Coordinate_Z_AP=[0];
%% SIM
Coordinate_X_SIM=0; 
Coordinate_Y_SIM=0;
Coordinate_Z_SIM=0;
SIM_Num=100; %
SIM_Thickness=5*lambda;
dx_mete_atoms=0.5*lambda;
dy_mete_atoms=0.5*lambda;
d_element_spacing = 0.5*lambda; %% Element spacing
SIM_layersV=[1:1:8];%8:-1:1
Sum_Maxmin_rate=0;
Sum_Sum_rate=0;
Sum_GM_rate=0;
%% Variable that need to be optimizaed
Power = 10.^((10)/10);%
%% Loop
for variable_n=1:1:length(SIM_layersV)
    Power_now=Power*10^(5/10);
    SIM_layers=SIM_layersV(variable_n);
    d_layer_spacing = SIM_Thickness/SIM_layers; %% Adjacent layer spacing in TX-SIM
    W_T_1 = zeros(SIM_Num,Antenna_Num);
    W_T = zeros(SIM_Num,SIM_Num);
    Corr_T = zeros(SIM_Num,SIM_Num);
    for SIM_n1 = 1:SIM_Num
        SIM_z1 = ceil(SIM_n1/sqrt(SIM_Num)); %% 
        SIM_x1 = mod(SIM_n1-1,sqrt(SIM_Num))+1; %% 
        % The channel from transmitter to the first layer of SIM
        for Antenna_n = 1:Antenna_Num
            d_transmit_2_SIM = sqrt(d_layer_spacing^2 + ...
                ( (SIM_x1-(1+sqrt(SIM_Num))/2)*d_element_spacing )^2 + ...
                ( (SIM_z1-(1+sqrt(SIM_Num))/2)*d_element_spacing - (Antenna_n-(1+Antenna_Num)/2)*d_antenna_spacing )^2 ); %% 
            %W_T_1(SIM_n1,Antenna_n) = lambda/4/pi/d_transmit_2_SIM*exp(-1i*2*pi*d_transmit_2_SIM/lambda);
            W_T_1(SIM_n1,Antenna_n) = lambda^2/4*(d_layer_spacing/d_transmit_2_SIM/d_transmit_2_SIM*(1/2/pi/d_transmit_2_SIM-1i/lambda))*exp(1i*2*pi*d_transmit_2_SIM/lambda);

        end
        % Calculate inter-layer transmission coefficient matrix W_T and channel correlation matrix Corr_T associated with TX-SIM
        for SIM_n2 = 1:SIM_Num
            SIM_z2 = ceil(SIM_n2/sqrt(SIM_Num)); %% 
            SIM_x2 = mod(SIM_n2-1,sqrt(SIM_Num))+1; %% 
            d_n1_2_n2_SameLayer= sqrt(  (SIM_x1-SIM_x2)^2 +  (SIM_z1-SIM_z2) ^2 )*d_element_spacing; %% 同一层元素n1到元素n2的距离
            d_n1_2_n2_AdjacentLayer = sqrt(d_layer_spacing^2 + d_n1_2_n2_SameLayer^2); %% 
            %W_T(SIM_n2,SIM_n1) = lambda/4/pi/d_n1_2_n2_AdjacentLayer*exp(-1i*2*pi*d_n1_2_n2_AdjacentLayer/lambda); %% old model
            W_T(SIM_n2,SIM_n1) = (dx_mete_atoms*dy_mete_atoms)*(d_layer_spacing/d_n1_2_n2_AdjacentLayer)/d_n1_2_n2_AdjacentLayer*(1/2/pi/d_n1_2_n2_AdjacentLayer-1i/lambda)*exp(1i*2*pi*d_n1_2_n2_AdjacentLayer/lambda); %% new model
            Corr_T(SIM_n2,SIM_n1) = sinc(2*d_n1_2_n2_SameLayer/lambda); %% 
        end
    end

    d_SIM_2_User=zeros(User_Num,1);
    PL_SIM_2_User=zeros(User_Num,1);
    alpha_RU=3.5; 
    for User_n=1:User_Num
        d_SIM_2_User(User_n)=sqrt( (Coordinate_X_User(User_n)-Coordinate_X_SIM)^2 +...
            (Coordinate_Y_User(User_n)-Coordinate_Y_SIM)^2 +...
            (Coordinate_Z_User(User_n)-Coordinate_Z_SIM)^2 );
        PL_SIM_2_User(User_n)=(lambda/(4*pi))^2*d_SIM_2_User(User_n)^(-alpha_RU);
    end
    
    %% MonteCarlo
    for monteCarlo_n=1:MonteCarlo
        % Channels from SIM to users
        Channel_SIM_User=zeros(SIM_Num,User_Num);
        for User_n=1:User_Num
            Channel_SIM_User(:,User_n)=sqrt(PL_SIM_2_User(User_n))*Corr_T^(1/2)*sqrt(1/2)*(randn(SIM_Num,1)+1i*randn(SIM_Num,1))*10^Enlarge;
        end
        % Initialize SIM phase shifts
        Theta=zeros(SIM_Num,SIM_layers);
        for SIM_layers_n = 1:SIM_layers
            omega = 0 + (2*pi-0)*rand(SIM_Num,1);
            Theta(:,SIM_layers_n)=(cos(omega)+1i*sin(omega));
        end

        SIM_Equivalent_Coefficient=diag(Theta(:,1));
        for SIM_layers_n = 2:SIM_layers
            SIM_Equivalent_Coefficient = diag(Theta(:,SIM_layers_n))*W_T*SIM_Equivalent_Coefficient;
        end

        Channel_Antenna_SIM_User=zeros(Antenna_Num,User_Num);
        for User_n=1:User_Num
            for Antenna_n = 1:Antenna_Num
                Channel_Antenna_SIM_User(Antenna_n,User_n)=Channel_SIM_User(:,User_n)'*SIM_Equivalent_Coefficient*W_T_1(:,Antenna_n);
            end
        end

        
        % Initialize power coefficient
        Power_EachAntenna=zeros(Antenna_Num,1);
        SqrtPower_EachAntenna=zeros(Antenna_Num,1);
        for Antenna_n = 1:Antenna_Num
            Power_EachAntenna(Antenna_n)=Power_now/Antenna_Num;
            SqrtPower_EachAntenna(Antenna_n)=sqrt(Power_EachAntenna(Antenna_n));
        end

        Maxmin_before=0;
        Maxmin_after=0;
        Middle_rate=0;
        loop_index=0;
        now_Convergence_coefficient=100;

        Error_xx=zeros(Antenna_Num,User_Num);Error_mu=zeros(1,User_Num);
        xx=zeros(Antenna_Num,User_Num);
        xx_2=zeros(Antenna_Num,User_Num);
        Q_UU=zeros(Antenna_Num,User_Num);
        Q_U=zeros(User_Num);
        Constant_c=zeros(User_Num,1);
        Penalty_Factor_BF=20;%2*sqrt(User_Num)
        allloop=1;allloop1=1;
        %% 
        while(abs(now_Convergence_coefficient)>=10^(-5))

            Channel_Multiply_BF=zeros(Antenna_Num,User_Num);
            Channel_Multiply_BF_2=zeros(Antenna_Num,User_Num);
            Channel_Multiply_Channel=zeros(Antenna_Num,User_Num);
            Rate=zeros(1,User_Num);
            Inverse_GM_Rate=zeros(1,User_Num);
            P=[1,0];
            Denominator=zeros(1,User_Num); 
            A=zeros(2,2,User_Num);
            B=zeros(2,2,User_Num);
            B11=zeros(1,User_Num);
            B12=zeros(1,User_Num);
            B21=zeros(1,User_Num);
            B22=zeros(1,User_Num);
            power_consume=zeros(1,Antenna_Num);

            for User_n=1:User_Num 
                for Antenna_n = 1:Antenna_Num %
                    Channel_Multiply_BF(Antenna_n,User_n)=Channel_Antenna_SIM_User(Antenna_n,User_n)*SqrtPower_EachAntenna(Antenna_n);
                    Channel_Multiply_BF_2(Antenna_n,User_n)=Channel_Multiply_BF(Antenna_n,User_n)*Channel_Multiply_BF(Antenna_n,User_n)';
                end
            end

            for User_n=1:User_Num
                Denominator(User_n)=sum(Channel_Multiply_BF_2(:,User_n))-Channel_Multiply_BF_2(User_n,User_n)+sigma;%速率表达式的分母项
                A(:,:,User_n)=[1,Channel_Multiply_BF(User_n,User_n)';
                    Channel_Multiply_BF(User_n,User_n),Denominator(User_n)+Channel_Multiply_BF_2(User_n,User_n)];
                B(:,:,User_n)=inv(A(:,:,User_n))*P'*inv(P*inv(A(:,:,User_n))*P')*P*inv(A(:,:,User_n));
                B11(User_n)=B(1,1,User_n);
                B12(User_n)=B(1,2,User_n);
                B21(User_n)=B(2,1,User_n);
                B22(User_n)=B(2,2,User_n);
                Constant_c(User_n)= real(log(det(P*inv(A(:,:,User_n))*P')) + trace(B(:,:,User_n)*A(:,:,User_n)) - ...
                    ( B11(User_n)+B22(User_n)*sigma ));
                Rate(User_n)=log(det(P*inv(A(:,:,User_n))*P'))/log(2);
                Rate2(User_n)=log(1+Channel_Multiply_BF_2(User_n,User_n)/Denominator(User_n))/log(2);
            end
            for User_n=1:User_Num
                Inverse_GM_Rate(User_n)=1/log(2);
            end

            sum(Rate)
            Maxmin_after= sum(Rate);
            all_loop_rate(allloop)=Maxmin_after;
                allloop=allloop+1;
            now_Convergence_coefficient=abs(Maxmin_after-Maxmin_before)/Maxmin_before;
            Maxmin_before=Maxmin_after;

            % Power_EachAntenna
            lanmuda_WW=1000*ones(1,1);
            lanmuda_up_WW=1*ones(1,1);
            lanmuda_down_WW=1*zeros(1,1);

            for User_n=1:User_Num 
                for Antenna_n = 1:Antenna_Num %
                    Q_U(Antenna_n)=-Inverse_GM_Rate(Antenna_n)*...
                        real(B12(Antenna_n)*Channel_Antenna_SIM_User(Antenna_n,Antenna_n));
                    Q_UU(Antenna_n,User_n)= real( Inverse_GM_Rate(User_n)*B22(User_n)*...
                        ( Channel_Antenna_SIM_User(Antenna_n,User_n)*Channel_Antenna_SIM_User(Antenna_n,User_n)' ) );
                end
            end
            for Antenna_n = 1:Antenna_Num %
                SqrtPower_EachAntenna(Antenna_n)=Q_U(Antenna_n)/sum(Q_UU(Antenna_n,:));
                power_consume(Antenna_n)=SqrtPower_EachAntenna(Antenna_n)^2;
            end
            if  sum(power_consume)>Power_now
                for Antenna_n = 1:Antenna_Num %
                    SqrtPower_EachAntenna(Antenna_n)=Q_U(Antenna_n)/(sum(Q_UU(Antenna_n,:))+lanmuda_WW*log(2));
                    power_consume(Antenna_n)=SqrtPower_EachAntenna(Antenna_n)^2;
                end
                if sum(power_consume)>Power_now
                    while (sum(power_consume)-Power_now>0)
                        lanmuda_WW=lanmuda_WW*100;
                        for Antenna_n = 1:Antenna_Num %
                            SqrtPower_EachAntenna(Antenna_n)=Q_U(Antenna_n)/(sum(Q_UU(Antenna_n,:))+lanmuda_WW*log(2));
                            power_consume(Antenna_n)=SqrtPower_EachAntenna(Antenna_n)^2;
                        end
                    end
                end
                lanmuda_up_WW=lanmuda_WW;
                lanmuda_down_WW=10^(-7);
                while (1)
                    lanmuda_WW=0.5*(lanmuda_up_WW+lanmuda_down_WW);
                    for Antenna_n = 1:Antenna_Num %
                        SqrtPower_EachAntenna(Antenna_n)=Q_U(Antenna_n)/(sum(Q_UU(Antenna_n,:))+lanmuda_WW*log(2));
                        power_consume(Antenna_n)=SqrtPower_EachAntenna(Antenna_n)^2;
                    end
                    if  sum(power_consume)>Power_now
                        lanmuda_down_WW = lanmuda_WW;
                    elseif  sum(power_consume)<Power_now
                        lanmuda_up_WW = lanmuda_WW ;
                    end
                    if abs(lanmuda_up_WW-lanmuda_down_WW)/lanmuda_down_WW<=10^(-5) ||  lanmuda_WW<10^(-7)
                        break;
                    end
                end
            end
            sum_power_consume=sum(power_consume);


            % Thet SIM_phase shift
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                X_SIM_Equivalent_Coefficient = zeros(SIM_Num,SIM_Num,SIM_layers);
                Y_SIM_Equivalent_Coefficient = zeros(SIM_Num,SIM_Num,SIM_layers);
                X_temp=eye(SIM_Num);
                Y_temp=eye(SIM_Num);
                for SIM_layers_n = 1:SIM_layers
                    if SIM_layers_n == 1
                        X_temp=eye(SIM_Num);
                        for SIM_layers_nn = 2 : SIM_layers
                            X_temp  =  diag(Theta(:,SIM_layers_nn))*W_T*X_temp;
                        end
                        X_SIM_Equivalent_Coefficient(:,:,SIM_layers_n) = X_temp;
                        Y_SIM_Equivalent_Coefficient(:,:,SIM_layers_n) = eye(SIM_Num);
                    elseif SIM_layers_n == SIM_layers
                        Y_temp=diag(Theta(:,1));
                        for SIM_layers_nn = 2 : SIM_layers-1
                            Y_temp  =  diag(Theta(:,SIM_layers_nn))*W_T*Y_temp;
                        end
                        Y_temp=W_T*Y_temp;
                        X_SIM_Equivalent_Coefficient(:,:,SIM_layers_n) = eye(SIM_Num);
                        Y_SIM_Equivalent_Coefficient(:,:,SIM_layers_n) = Y_temp;
                    else
                        X_temp=eye(SIM_Num);
                        for SIM_layers_nn = SIM_layers_n+1 : SIM_layers
                            X_temp  =  diag(Theta(:,SIM_layers_nn))*W_T*X_temp;
                        end
                        X_SIM_Equivalent_Coefficient(:,:,SIM_layers_n) = X_temp;

                        Y_temp=diag(Theta(:,1));
                        for SIM_layers_nn = 2 : SIM_layers_n-1
                            Y_temp  =  diag(Theta(:,SIM_layers_nn))*W_T*Y_temp;
                        end
                        Y_temp=W_T*Y_temp;
                        Y_SIM_Equivalent_Coefficient(:,:,SIM_layers_n) = Y_temp;
                    end
                    %aaaaaa=    X_SIM_Equivalent_Coefficient(:,:,SIM_layers_n) * diag(Theta(:,SIM_layers_n)) * Y_SIM_Equivalent_Coefficient(:,:,SIM_layers_n);

                    D_temp1=zeros(SIM_Num,SIM_Num,Antenna_Num,User_Num);
                    V_temp1=zeros(SIM_Num,1,User_Num);
                    for User_n=1:User_Num %第几个用户
                        m_temp=Inverse_GM_Rate(User_n)*B12(User_n)*Channel_SIM_User(:,User_n)'*X_SIM_Equivalent_Coefficient(:,:,SIM_layers_n);
                        n_temp=Y_SIM_Equivalent_Coefficient(:,:,SIM_layers_n)*W_T_1(:,User_n)*SqrtPower_EachAntenna(User_n);
                        V_temp1(:,:,User_n)=(m_temp.*n_temp.');
                        for Antenna_n = 1:Antenna_Num %
                            mm_temp=Channel_SIM_User(:,User_n)'*X_SIM_Equivalent_Coefficient(:,:,SIM_layers_n);
                            nn_temp=Y_SIM_Equivalent_Coefficient(:,:,SIM_layers_n)*W_T_1(:,Antenna_n)*SqrtPower_EachAntenna(Antenna_n);   
                            D_temp1(:,:,Antenna_n,User_n)=Inverse_GM_Rate(User_n)*B22(User_n)*(diag(mm_temp)*nn_temp*nn_temp'*diag(mm_temp)').';
                        end
                    end

                    D_temp=sum(sum(D_temp1,4),3);
                    V_temp=(sum(V_temp1,3));

                    Part_1=zeros(SIM_Num,1);
                    Part_2=zeros(SIM_Num,1);
                    for n=1:SIM_Num
                        Part_1(n)=V_temp(n);
                        for m=1:SIM_Num
                            if m==n
                                Part_2(n)=Part_2(n)+0;
                            else
                                Part_2(n) = Part_2(n) + D_temp( m , n ) * conj(Theta(m,SIM_layers_n));
                            end
                        end
                        Theta(n,SIM_layers_n)=exp(1i*(pi-angle(Part_1(n)+Part_2(n))));
                    end
                end

                SIM_Equivalent_Coefficient=diag(Theta(:,1));
                for SIM_layers_n = 2:SIM_layers
                    SIM_Equivalent_Coefficient = diag(Theta(:,SIM_layers_n))*W_T*SIM_Equivalent_Coefficient;
                end
                Channel_Antenna_SIM_User=zeros(Antenna_Num,User_Num);%横坐标是用户索引
                for User_n=1:User_Num
                    for Antenna_n = 1:Antenna_Num
                        Channel_Antenna_SIM_User(Antenna_n,User_n)=Channel_SIM_User(:,User_n)'*SIM_Equivalent_Coefficient*W_T_1(:,Antenna_n);
                    end
                end
            loop_index=loop_index+1;
            ;
        end
        disp(['Next Iteration: ']);
        Sum_Maxmin_rate=Sum_Maxmin_rate+min(Rate);
        Sum_Sum_rate=Sum_Sum_rate+sum(Rate);
        Sum_GM_rate=Sum_GM_rate+Maxmin_after;
        ;
    end
    Ave_Maxmin_rate(variable_n)=Sum_Maxmin_rate/(MonteCarlo);
    Ave_Sum_rate(variable_n)=Sum_Sum_rate/(MonteCarlo);
    Ave_GM_rate(variable_n)=Sum_GM_rate/(MonteCarlo);
    Sum_Maxmin_rate=0;
    Sum_Sum_rate=0;
    Sum_GM_rate=0;
    Maxmin_after=0;
end
code_time_secend=toc;
code_time_min=toc/60;