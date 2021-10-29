classdef TFM_ANALYSIS < handle
    properties (SetAccess = private, Hidden = true)
        parameters;
        utility;
    end
    methods(Static)
        function params=get_default_parameters(init_params)
            %OPTICAL PARAMETERS
%             params=BASIC_OPTICAL_PARAMETER();
            %SIMULATION PARAMETERS
            params.resolution=[1 1 1];
            params.use_GPU = true;
            % Subset spacing (pix) 
            w0 = 32;
            d0 = 6;
            params.blocksizes = [w0 w0 8]; % [iblocksize jblocksize,kblocksize]
            params.overlap = 1-d0/w0;
            params.padding = 16;
            params.method = '3Dcentroid';
            params.N = 2;
            params.Ni = 4;
            
            if nargin==1
                params=update_struct(params,init_params);
            end
        end
    end
    methods
        function h=TFM_ANALYSIS(params)
            h.parameters=params;
        end

        function Quiver = PIV_3D(h,im1,im2)

        % 1. GPU Initialization
            if h.parameters.use_GPU
                im1 = gpuArray(single(im1));
                im2 = gpuArray(single(im2));
            end
            h.parameters.sizes = size(im1, 1:3);
            if any(size(im1) ~= size(im2))
                error('Image dimensions do not match together.')
            end
            
        % 2. Padarray; xExpandMatrix is equivalent to padarray.
            if size(im1,3) > 1
                im1 = padarray(im1, h.parameters.padding*ones(1,3), 0);
                im2 = padarray(im2, h.parameters.padding*ones(1,3), 0);
            else
                im1 = padarray(im1, h.parameters.padding*ones(1,2), 0);
                im2 = padarray(im2, h.parameters.padding*ones(1,2), 0);
            end
            %{
            im1 = xExpandMatrix(im1, 1, 1, 1, padding, padding, padding, padding, padding, padding); % pad
            im2 = xExpandMatrix(im2, 1, 1, 1, padding, padding, padding, padding, padding, padding);
            %}

        % 3. Define tile increments
            incs = round(h.parameters.blocksizes* (1-h.parameters.overlap)); %[inci, incj, inck] - define increment
            if sum((incs < 1) + (incs > h.parameters.blocksizes)) > 0 
                error('Wrong Overlap in Correlation Algorithm')
            end
            %{
            inci = round(iblocksize*(1-overlap)); % define increment
            incj = round(jblocksize*(1-overlap)); % define increment
            inck = round(kblocksize*(1-overlap)); % define increment
            %}
            
        % 4. Initialize Quiver
            Quiver = struct;
            if size(im1,3) > 1
                Qsz = floor((h.parameters.sizes - h.parameters.blocksizes+1-1) ./ incs)+1;
            else
                Qsz = floor((h.parameters.sizes(1:2) - h.parameters.blocksizes(1:2)+1-1) ./ incs(1:2))+1;
                Qsz(3) = 1;
            end
            Quiver.U = zeros(Qsz,'single');
            Quiver.V = zeros(Qsz,'single');
            Quiver.W = zeros(Qsz,'single');
            if length(h.parameters.sizes(3)) == 2
                Quiver.Z = zeros(Qsz,'single')+inf;
                Quiver.W = zeros(Qsz,'single')+inf;
            end
            Quiver.pkh = zeros(Qsz,'single');

            ki = 1; kj = 1; kk = 1;

            if Qsz(3) > 1
                im11 = im1(ki*(incs(1)-1)+1:ki*(incs(1)-1)+h.parameters.blocksizes(1)+2*h.parameters.padding,...
                    kj*(incs(2)-1)+1:kj*(incs(2)-1)+h.parameters.blocksizes(2)+2*h.parameters.padding,...
                    kk*(incs(3)-1)+1:kk*(incs(3)-1)+h.parameters.blocksizes(3)+2*h.parameters.padding);
                im11 = im11(h.parameters.padding+1:end-h.parameters.padding,...
                    h.parameters.padding+1:end-h.parameters.padding,...
                    h.parameters.padding+1:end-h.parameters.padding);
            else
                im11 = im1(ki*(incs(1)-1)+1:ki*(incs(1)-1)+h.parameters.blocksizes(1)+2*h.parameters.padding,...
                    kj*(incs(2)-1)+1:kj*(incs(2)-1)+h.parameters.blocksizes(2)+2*h.parameters.padding,:);
                im11 = im11(h.parameters.padding+1:end-h.parameters.padding,...
                    h.parameters.padding+1:end-h.parameters.padding,:);
            end
            [h.utility.YY,h.utility.XX, h.utility.ZZ] = ndgrid(1:size(im11,1), 1:size(im11,2),1:size(im11,3));
            if h.parameters.use_GPU
                h.utility.XX = gpuArray(h.utility.XX);
                h.utility.YY = gpuArray(h.utility.YY);
                h.utility.ZZ = gpuArray(h.utility.ZZ);
            end

            %{
            Size1 = size([1:inci:sz(1)-iblocksize+1],2); % DT
            Size2 = size([1:incj:sz(2)-jblocksize+1],2); % DT
            Size3 = size([1:inck:sz(3)-kblocksize+1],2); 
            x =   zeros(Size1,Size2,Size3); % position in x
            y =   zeros(Size1,Size2,Size3); % position in y
            z =   zeros(Size1,Size2,Size3); % position in z
            dx =  zeros(Size1,Size2,Size3); % displacements in x
            dy =  zeros(Size1,Size2,Size3); % displacements in y
            dz =  zeros(Size1,Size2,Size3); % displacements in z
            pkh = zeros(Size1,Size2,Size3); % height of the xcorr peak
            %}

        % 5. Major Loop
            tic;           
            for ki = 1:Qsz(1)
                for kj = 1:Qsz(2)
                    for kk = 1:Qsz(3)
                        clc,
                        disp([num2str(ki) ' / ' num2str(Qsz(1))])
                        disp([num2str(kj) ' / ' num2str(Qsz(2))])
                        disp([num2str(kk) ' / ' num2str(Qsz(3))])
            
                        niterations=0; DX=inf; DY=inf; DZ=inf; % initialize iterative process
                        while niterations<=h.parameters.N && abs(Quiver.U(ki, kj, kk)-DX) > 0.02 &&...
                                abs(Quiver.V(ki, kj, kk)-DY) > 0.02 && abs(Quiver.W(ki, kj, kk)-DZ) > 0.02
                        
                        % 5-1. allocate the displacement on the new matrix
                            niterations = niterations + 1;
                            DX = Quiver.U(ki, kj, kk);
                            DY = Quiver.V(ki, kj, kk);
                            DZ = Quiver.W(ki, kj, kk);
                            %{
                                DX=(dx((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck));
                                DY=(dy((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck));
                                DZ=(dz((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck)); 
                            %}
            
                        % 5-2. Initialize block
                            if ~(DX == 0 && DY == 0) %skip first iteration(DX=0,Dy=0 for first iteration)
                                im23 = circshift_subpixel_MS(im22,[DY DX DZ],h.parameters.use_GPU);
                            else
                                if Qsz(3) > 1
                                    im11 = im1(ki*(incs(1)-1)+1:ki*(incs(1)-1)+h.parameters.blocksizes(1)+2*h.parameters.padding,...
                                        kj*(incs(2)-1)+1:kj*(incs(2)-1)+h.parameters.blocksizes(2)+2*h.parameters.padding,...
                                        kk*(incs(3)-1)+1:kk*(incs(3)-1)+h.parameters.blocksizes(3)+2*h.parameters.padding);
                                    im22 = im2(ki*(incs(1)-1)+1:ki*(incs(1)-1)+h.parameters.blocksizes(1)+2*h.parameters.padding,...
                                        kj*(incs(2)-1)+1:kj*(incs(2)-1)+h.parameters.blocksizes(2)+2*h.parameters.padding,...
                                        kk*(incs(3)-1)+1:kk*(incs(3)-1)+h.parameters.blocksizes(3)+2*h.parameters.padding);
                                    im11 = im11(h.parameters.padding+1:end-h.parameters.padding,...
                                        h.parameters.padding+1:end-h.parameters.padding,...
                                        h.parameters.padding+1:end-h.parameters.padding);


                                else
                                    im11 = im1(ki*(incs(1)-1)+1:ki*(incs(1)-1)+h.parameters.blocksizes(1)+2*h.parameters.padding,...
                                        kj*(incs(2)-1)+1:kj*(incs(2)-1)+h.parameters.blocksizes(2)+2*h.parameters.padding,:);
                                    im22 = im2(ki*(incs(1)-1)+1:ki*(incs(1)-1)+h.parameters.blocksizes(1)+2*h.parameters.padding,...
                                        kj*(incs(2)-1)+1:kj*(incs(2)-1)+h.parameters.blocksizes(2)+2*h.parameters.padding,:);
                                    im11 = im11(h.parameters.padding+1:end-h.parameters.padding,...
                                        h.parameters.padding+1:end-h.parameters.padding,:);
                                end
                                im23 = im22;
                            end
                            if Qsz(3) > 1
                                im23 = im23(h.parameters.padding+1:end-h.parameters.padding,...
                                    h.parameters.padding+1:end-h.parameters.padding,...
                                    h.parameters.padding+1:end-h.parameters.padding);
                            else
                                im23 = im23(h.parameters.padding+1:end-h.parameters.padding,...
                                    h.parameters.padding+1:end-h.parameters.padding,:);
                            end
                            %{
                                im11 = im1(ki : ki+iblocksize+2*padding-1 , kj : kj+jblocksize+2*padding-1, kk:kk+kblocksize+2*padding-1); % crop the block with padding
                                im22 = im2(ki : ki+iblocksize+2*padding-1 , kj : kj+jblocksize+2*padding-1, kk:kk+kblocksize+2*padding-1);
                                im11=im11/mean2(im11);
                                im22=im22/mean2(im22);
                                im11 = xsubpix_shift(im11,DX/2,DY/2,DZ/2); % subpixel shift of the image
                                im22 = xsubpix_shift(im22,-DX/2,-DY/2,-DZ/2);
                                im11 = im11(padding+1 : padding+iblocksize , padding+1 : padding+jblocksize, padding+1 : padding+kblocksize); % crop the block
                                im22 = im22(padding+1 : padding+iblocksize , padding+1 : padding+jblocksize, padding+1 : padding+kblocksize);
                                im11 = im11 - mean2(im11); % subtract the mean
                                im22 = im22 - mean2(im22);    
                          % I won't use window function
                                if(window),   
                                    im11 = xwindow2D(im11,window); % multiply each block by a window
                                    im22 = xwindow2D(im22,window);                
                                end
                            %}

                    % 5-3. Get correlation & quiver vectors
                        [~, dshift, Quiver.pkh(ki,kj,kk)] = h.Get_correlation_3D(im11, im23,h.parameters.N);
                        %%
                        Quiver.V(ki, kj, kk) = Quiver.V(ki, kj, kk) + dshift(2);
                        Quiver.U(ki, kj, kk) = Quiver.U(ki, kj, kk) + dshift(1);
                        Quiver.W(ki, kj, kk) = Quiver.W(ki, kj, kk) + dshift(3);
                        %{
                        c = xcorrf2(im22,im11,'yes'); % / (std(im11(:))*std(im22(:))*blocksize^2); % compute the correlation function in Fourier Space and normalize
                        c = real(c(1:end-1,1:end-1,1:end-1)); % resize for x_cntr 
                        x((ki+inci-1)/inci, (kj+incj-1)/incj , (kk+inck-1)/inck) = kj+jblocksize/2;
                        y((ki+inci-1)/inci, (kj+incj-1)/incj , (kk+inck-1)/inck) = ki+iblocksize/2;
                        z((ki+inci-1)/inci, (kj+incj-1)/incj , (kk+inck-1)/inck) = kk+kblocksize/2;
                        dx((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck) = DX + x_cntr(c,N,method,'X',th); % compute the displacement
                        dy((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck)  = DY + x_cntr(c,N,method,'Y',th);   
                        dz((ki+inci-1)/inci, (kj+incj-1)/incj,  (kk+inck-1)/inck) = DZ + x_cntr(c,N,method,'Z',th); 
                        pkh((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck) = max3(c); % store peak height
                        %}
            
                        end
                    end
                end
            end
            Quiver.U = gather(Quiver.U);
            Quiver.V = gather(Quiver.V);
            Quiver.W = gather(Quiver.W);

            if length(Qsz) > 2
            [Quiver.Y Quiver.X Quiver.Z] = ndgrid((1:Qsz(1))*(h.parameters.blocksizes(1)-1)+h.parameters.blocksizes(1)/2,...
                (1:Qsz(2))*(h.parameters.blocksizes(2)-1)+h.parameters.blocksizes(2)/2,...
                (1:Qsz(3))*(h.parameters.blocksizes(3)-1)+h.parameters.blocksizes(3)/2);
            else
            [Quiver.Y Quiver.X Quiver.Z] = ndgrid((1:Qsz(1))*(h.parameters.blocksizes(1)-1)+h.parameters.blocksizes(1)/2,...
                (1:Qsz(2))*(h.parameters.blocksizes(2)-1)+h.parameters.blocksizes(2)/2,1);
            end
            Quiver.tm = toc;
        end

        function [B] = circshift_subpixel_MS(h,B,shift0)
%             [YY, XX, ZZ] = ndgrid(1:size(B,1), 1:size(B,2),1:size(B,3));
%             if Bool_GPU
%                 XX = gpuArray(XX);
%                 YY = gpuArray(YY);
%                 ZZ = gpuArray(ZZ);
%             end
            if all(mod(shift0,1) == 0) % integer shift
                B = circshift(B, shift0);
            else % subpixel shift
                B =  fftshift(fftn(ifftshift(B))).* ...
                    exp(-1i.*2*pi*(shift0(1).*h.utility.YY/size(B,1)+shift0(2).*h.utility.XX/size(B,2)+shift0(3).*h.utility.ZZ/size(B,3)));
                B = fftshift(ifftn(ifftshift(B)));
                B = abs(B);
            end
        end

         function [Corr3D, R, cpeak] = Get_correlation_3D(h,im11, im22,N)
            im11 = (im11-mean(im11(:))) ./ std(im11(:),0);
            im22 = (im22-mean(im22(:))) ./ std(im22(:),0);
        
            Corr3D = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(im11))).*...
                conj(fftshift(fftn(ifftshift(im22)))))))) / prod(size(im11));
            R = h.peak_subpixel_positioner(Corr3D,N);
            cpeak = max(Corr3D(:));
         end       

        function R = peak_subpixel_positioner(h,Corr3D,N)
           [mxv,idx] = max(Corr3D(:));
           [mi,mj,mk] = ind2sub(size(Corr3D),idx);
           if length(size(Corr3D)) > 2
                       % image near the maximum value
                    if mi-N-1>0 && mi+N+1<size(Corr3D,1) && mj-N-1>0 && mj+N+1<size(Corr3D,2) ...
                            && mk-N-1>0 && mk+N+1<size(Corr3D,3)
                    
                        im_1 = Corr3D(mi-N-1:mi+N+1,mj-N-1:mj+N+1,mk-N-1:mk+N+1);
                        % Matrix without maximum region
                        im_wm = Corr3D;
                        im_wm(mi-N-1:mi+N+1,mj-N-1:mj+N+1,mk-N-1:mk+N+1) = 0;
                        th = max(im_wm(:));
                        im_1 = im_1-th;
                        im_1(im_1<0)=0;
                        
                        meshi = mi-N-1 : mi+N+1;
                        meshj = mj-N-1 : mj+N+1;
                        meshk = mk-N-1 : mk+N+1;
                        
                        [grdj,grdi,grdk] = ndgrid(meshi,meshj,meshk); 
                        
                        xc=sum(im_1(:).*grdi(:))/sum(im_1(:));
                        yc=sum(im_1(:).*grdj(:))/sum(im_1(:));
                        zc=sum(im_1(:).*grdk(:))/sum(im_1(:));
                    else
                         xc=size(Corr3D,2)/2-1;    
                         yc=size(Corr3D,1)/2-1;    
                         zc=size(Corr3D,3)/2-1;
                    end
           else
                    if mi-N-1>0 && mi+N+1<size(Corr3D,1) && mj-N-1>0 && mj+N+1<size(Corr3D,2)
                    
                        im_1 = Corr3D(mi-N-1:mi+N+1,mj-N-1:mj+N+1);
                        im_wm = Corr3D;
                        im_wm(mi-N-1:mi+N+1,mj-N-1:mj+N+1) = 0;
                        th = max(im_wm(:));
                        % Matrix without maximum region
                        im_1(im_1<0)=0;
                        
                        meshi = mi-N-1 : mi+N+1;
                        meshj = mj-N-1 : mj+N+1;
                        
                        [grdj,grdi] = ndgrid(meshi,meshj); 
                        
                        xc=sum(im_1(:).*grdi(:))/sum(im_1(:));
                        yc=sum(im_1(:).*grdj(:))/sum(im_1(:));
                        zc=1;
                    else
                         xc=size(Corr3D,2)/2-1;    
                         yc=size(Corr3D,1)/2-1;    
                         zc=size(Corr3D,3)/2-1;
                    end
           end
            R = [xc yc zc] - floor(size(Corr3D,[2 1 3])/2) - 1;
            R = R .* h.parameters.resolution;
        %     R
        %      figure(1),imagesc(max(Corr3D,[],3)), colorbar
        end

        
        function [H,XX,YY,ZZ] = mk_ellipse_MS(h,sizes, radii)
            if length(radii) == 1
                radii = ones(1,3) * radii;
            end
            if sizes(3) > 1
                [YY,XX,ZZ] = ndgrid(1:sizes(1), 1:sizes(2), 1:sizes(3));
                XX = XX - floor(sizes(2)/2) - 1;
                YY = YY - floor(sizes(1)/2) - 1;
                ZZ = ZZ - floor(sizes(3)/2) - 1;
                H = (XX./radii(2)).^2+(YY./radii(1)).^2+(ZZ./radii(3)).^2<=1.0;
            else
                ZZ = [];
                [YY,XX] = ndgrid(1:sizes(1), 1:sizes(2));
                XX = XX - floor(sizes(2)/2) - 1;
                YY = YY - floor(sizes(1)/2) - 1;
                H = (XX./radii(2)).^2+(YY./radii(1)).^2<=1.0;
            end
        end
            
    end

end


