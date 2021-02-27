function iter2dat(varargin)
% ITER2DAT A format converter which convert *.iter file and related mesh
% files to so called 'x,y,z' *.data file for post-processing.
%
% The *.iter file reading functions come from 'plotOccam2DMT.m' routine
% of SIO, UCSD.
%
%   ------------------
%   Bo Yang, 2011.
%   China Univ. of Geosciences, Wuhan, China.
%   Comments, bug reports and questions, please send to:
%   yangbo.cug@163.com.
%   Copyright 2011-2016 Bo Yang.
%   $Revision: 1.0 $ $Date: 2012/04/05 21:50:20 $

%   Revision log:
%   2012/04/04 : Version 1.0 released.
%
    x0 = 0;
    plot_zero       = 0;
    
    cFileSpec = {'*.iter','Iteration files (*.iter)';'*.*','All Files'};
    [sIterFile,sIterPath] ...
        = uigetfile( cFileSpec, 'Select OCCAM iteration file to plot:' );
    if ~ischar(sIterFile)
        return
    end
%--------------------------------------------------------------------------
% READ ITERATION FILE:
%--------------------------------------------------------------------------
    [stIter, nParamValues] = readITER( fullfile( sIterPath, sIterFile ) );
    
%--------------------------------------------------------------------------
% READ MODEL FILE:
%--------------------------------------------------------------------------
    [fidmod,modfrm,modnam,moddes,mshfile,mshtyp,stcfile,...
     prjfile,bndofs,nlayer,irz,nrcol,columns,nexep,l1,...
     c1,l2,c2,pen]= readMODEL( fullfile( sIterPath, stIter.sModelFile ) ); 
 
%--------------------------------------------------------------------------
% READ MESH FILE:
%--------------------------------------------------------------------------
    [fidmsh,mshdesc,idx,nodey,nodez,nres,nfre,nexec,...
     fixres,freqs,meshy,meshz,nrc,nrl,mshres]...
     = readMESH( fullfile( sIterPath, mshfile ) );

%--------------------------------------------------------------------------
% READ DATA FILE:
%--------------------------------------------------------------------------
    [datfrm,datdes,nsites,sites,offsts,nfre,freqs, ndblks,DATA] ...
        = readDATA( fullfile( sIterPath, stIter.sDataFile ) );
    
%--------------------------------------------------------------------------
% MAKE FE MESH ELEMENTS:
%--------------------------------------------------------------------------
    sumy = cumsum([0; meshy]);
    sumz = cumsum([0; meshz]);
    
    % ADJUST SUMY WITH BINDING OFFSET
    % BNDOFS is the right terminating edge of the left-most reg. brick
    temp = bndofs - sumy(columns(1,1)+1);
    sumy = sumy+temp;
    
    % remove x0 offset for plotting:
    sumy = sumy-x0 - plot_zero;
    sumy = sumy/1000;
    sumz = sumz/1000; % convert to km
    offsts = (offsts - plot_zero) / 1000;
    x0 = x0/1000;
    
    [YY, ZZ, restri] = makeELEMENTS(sumy,sort(-1*sumz),mshres);
    ZZ = -1*ZZ;
    
%--------------------------------------------------------------------------
% MAKE FE MESH LINES:
%--------------------------------------------------------------------------
	[YLy, YLz, ZLy, ZLz] = meshlines(sumy,sumz);
 % h_mesh=plot(YLy,YLz,'k-',ZLy,ZLz,'k-'); axis ij
 % to turn mesh off: set(h_mesh,'visible','off');
 % to trun on mesh:  set(h_mesh,'visible','on');
 
%--------------------------------------------------------------------------
% MAKE REGULARIZATION GRID
%--------------------------------------------------------------------------
    [REGY, REGZ, freebrick, seabrick, fixedbrick] ...
        = assembleREG(nlayer,nrcol,columns,irz,mshres,sumy,sumz); 

%--------------------------------------------------------------------------
% EXTRACT THE AREA, ADDED BY YANG BO
%--------------------------------------------------------------------------
    nmod = length(nParamValues);
    for k = 1:nmod
        for j = 1:4
            outmod((k-1)*4+j,1) = REGY(j,k);
            outmod((k-1)*4+j,2) = REGZ(j,k);
            outmod((k-1)*4+j,3) = nParamValues(k);
        end
    end
 
    dlg_title = 'Input the extract area';
    prompt = {'Ymin(km):','Ymax(km):','Zmin(km):','Zmax(km):'};
	num_lines= 1;
	def     = {'-10','520','0','200'};
	answer  = inputdlg(prompt,dlg_title,num_lines,def);
    ymin = str2num(answer{1});
    ymax = str2num(answer{2});
    zmin = str2num(answer{3});
    zmax = str2num(answer{4});
    k = 0;
    for j = 1:nmod*4
        yy = outmod(j,1);
        zz = outmod(j,2);
        cc = outmod(j,3);
        if (yy > ymin && yy < ymax && zz > zmin && zz < zmax)
            k = k + 1;
            dd(k,:) = [yy -zz cc];
        end
    end
    [filename,pathname] = uiputfile('.dat','Save XYZ file');
    pathfile = strcat(pathname,filename);
    save(pathfile,'dd','-ascii');

%--------------------------------------------------------------------------
% EXTRACT THE SITES INFO, ADDED BY YANG BO
%--------------------------------------------------------------------------
    fid = fopen([pathfile,'_sites.bln'],'w');
    % Find nearest node:
    for i = 1:nsites
        [temp,I]  = min(abs(sumy+x0-offsts(i)));
        % I = find(temp==min(temp));       
        ysites(i) = sumy(I);
%        disp( sprintf('%s %8g %8g %i',sites{i}, offsts(i), ysites(i), I) );
        J1 = max(find(mshres(1:4:end,I-1)=='Z'));  %top triangle
        J4 = max(find(mshres(4:4:end,I-1)=='Z'));  %right triangle
        if isempty(J4) % fix for apple mesh right side TP sites
            J4 = 0;
        end
        if isempty(J1)
            J1 = 0;
        end
        if J1>J4       %site is on top of element
             zsites(i) = (sumz(J1));              
        elseif J1==J4  %site on bottom of this element
             zsites(i) = sumz(J1+1);

        elseif J4>J1
            disp('Error in locating site on seafloor');
            disp('There seems to be a seawater cavity at y =')
            disp(offsts(i))
        end
%         h_sites(i) = text(ysites(i),zsites(i),'\diamondsuit',...
%                      'HorizontalAlignment','center',...
%                      'VerticalAlignment','baseline' );
%         if isempty( plt_axes )
%             set( h_sites(i), 'FontSize',20);
%         end
        fprintf(fid,'%f,%f,%s\n',ysites(i),zsites(i),sites{i});
    end
    fclose(fid);
%     if strcmpi(plt_sitescenter,'on') % plot near surface, model center
%         ax = axis;
%         rng = abs(max(ysites)-min(ysites));
%         pct = rng * .2;    % show an extra % on each side
%         axis( [min(ysites)-pct max(ysites)+pct 0 max(20,rng)]); 
%     end
    % PLOT SITE NAMES
%     if isempty( plt_axes )
%         h_sitesnames = text( ysites, zeros(size(zsites)), sites ...
%             , 'VerticalAlignment', 'top' ...
%             , 'HorizontalAlignment', 'center' ...
%             );
%     end


return



%--------------------------------------------------------------------------    
function [stIter, nParams] = readITER( sFile )
%--------------------------------------------------------------------------
% readITER - read & parse an Occam2D iteration file.
%
% v1.0, July 2001 by: Kerry Key IGPP/SIO/UCSD
% v2.0, Nov 2006 by David Myer - rewritten to support new OCCAMITER_FLEX
% format files which have new entries for diagonal penalties, and roughness
% modification.  And which no longer require entries to be in a particular
% order.
%
% Parameters:
%   sFile   - path+name of iteration file to open.
% Returns:
%   stIter  - structure with members for each of the header items.  See code
%           for a list of possible members.  Structure will only have members
%           for the items that exist in the iteration file.  So not all
%           members will be present every time.
%   nParams - parameter values
%--------------------------------------------------------------------------

    % Open the iteration file
    fid = fopen( sFile, 'r' );
    
    % Default a value or two that should always exist.
    stIter.nMisfit      = 1000.0;
    stIter.nRoughness   = 0.0;
    
    % The header of the file consists of a bunch of "code: value" pairs
    % terminated by a pair that tells how many parameters there are (and is
    % immediately followed by the parameter list).
    % Read & process the pairs.
    while ~ feof(fid)
        % Get the current line & break up the code / value pair
        sLine = fgets( fid );
        [sCode, sValue] = strtok( sLine, ':' );
        sCode = lower(strtrim(sCode));
        sValue(1) = [];     % remove the leading token
        
        % If there is a user comment in the value, eliminate it.
        sValue = strtrim( strtok(sValue, '!%') );
        
        % Which code do we have?
        switch (sCode)
        case {'format'}
            stIter.sFormat = sValue;
            
        case {'description'}
            stIter.sDesc = sValue;
            
        case {'model file'}
            stIter.sModelFile = sValue;
            
        case {'data file'}
            stIter.sDataFile = sValue;
            
        case {'date/time'}
            stIter.sDateTime = sValue;
            
        case {'max iter', 'iterations to run'}
            stIter.nIterToRun = str2double( sValue );
            
        case {'req tol', 'target misfit'}
            stIter.nTargetMisfit = str2double( sValue );
            
        case {'iruf', 'roughness type'}
            stIter.nRoughnessType = str2double( sValue );
            stIter.bDiagPenalties = (stIter.nRoughnessType < 0);
            stIter.nRoughnessType = abs( stIter.nRoughnessType );
            
        case {'diagonal penalties'}
            stIter.bDiagPenalties = (str2double( sValue ) ~= 0);
            if isfield( stIter, 'nRoughnessType' )
                stIter.nRoughnessType = abs( stIter.nRoughnessType );
            end
            
        case {'debug level'}
            stIter.DebugLevel = str2double( sValue );
            
        case {'iteration'}
            stIter.nIteration = str2double( sValue );
            
        case {'pmu', 'lagrange value'}
            stIter.nLagrange = str2double( sValue );
            
        case {'rlast', 'roughness value'}
            stIter.nRoughness = str2double( sValue );
            
        case {'tlast', 'misfit value'}
            stIter.nMisfit = str2double( sValue );
            
        case {'ifftol', 'misfit reached'}
            stIter.bReachedMisfit = (str2double( sValue ) ~= 0);
            
        case {'modify roughness'}
            stIter.sModifyRoughness = sValue;
            
        case {'no. parms', 'param count'}
            stIter.nParams = str2double( sValue );
            break;	% The next thing in the file is the data param list!
            
            otherwise % do nothing
                
%         case default
%             disp( 'Error reading startup/iteration file!' );
%             error( ['Unknown or unsupported code: ' sCode] );
            
        end
        
    end
    if feof(fid) 
        error( ['Error: iteration file ' sFile ' terminated early!'] );
    end
    
    % Read the parameters
    nParams = fscanf( fid, '%f', stIter.nParams );
    
    % Close up & return
    fclose(fid);

    return




%--------------------------------------------------------------------------
function [fidmod,modfrm,modnam,moddes,mshfile,mshtyp,stcfile,...
         prjfile,bndofs,nlayer,irz,nrcol,columns,nexep,l1,...
         c1,l2,c2,pen]= readMODEL(modfile)
%--------------------------------------------------------------------------
% function readMODEL.m
% Version 1.0, July 2001
% Written by: Kerry Key IGPP/SIO/UCSD
% 
% A matlab script to read an OCCAM2DMT MODEL file
%
% Usage: 
%        [fidmod,modfrm,modnam,moddes,mshfile,mshtyp,stcfile,...
%         prjfile,bndofs,nlayer,irz,nrcol,columns,nexep,l1,...
%         c1,l2,c2,pen]= readMODEL(modfile);
% 
% Inputs
%        modfile: OCCAM2DMT MODEL file name (string format)
%
% Outputs: 
% 		fidmod
% 		modfrm
% 		modnam
% 		moddes
% 		mshfile
% 		mshtyp
% 		stcfile
% 		prjfile
% 		bndofs
% 		nlayer
% 		irz(nlayer)
% 		nrcol(nlayer)
% 		columns(nlayer,nrcol(nlayer))
% 		nexep
% 		l1(nexep)
% 		c1(nexep)
% 		l2(nexep)
% 		c2(nexep)
% 		pen(nexep)
%
%--------------------------------------------------------------------------


% OPEN AND READ IN MODEL FILE:
 
   fidmod  = fopen(modfile,'r');
    cline  = fgets(fidmod);
   modfrm  = sscanf(cline,'FORMAT:           %80c');
    cline  = fgets(fidmod);
   modnam  = sscanf(cline,'MODEL NAME:       %80c');
    cline  = fgets(fidmod);
   moddes  = sscanf(cline,'DESCRIPTION:      %80c');
    cline  = fgets(fidmod);
   mshfile = sscanf(cline,'MESH FILE:        %s');
    cline  = fgets(fidmod);  
   mshtyp  = sscanf(cline,'MESH TYPE:        %s');
    cline  = fgets(fidmod); 
   stcfile = sscanf(cline,'STATICS FILE:     %s');
    cline  = fgets(fidmod);
   prjfile = sscanf(cline,'PREJUDICE FILE:   %s');
    cline  = fgets(fidmod);
   bndofs  = sscanf(cline,'BINDING OFFSET:   %f');
    cline  = fgets(fidmod);
   nlayer  = sscanf(cline,'NUM LAYERS:       %f');
   for i = 1:nlayer
       irz(i)   = fscanf(fidmod,'%i',1);
       nrcol(i) = fscanf(fidmod,'%i',1);
       if i==1;   % initialize columns array  
           columns = nan*ones(nlayer,nrcol(1));
       end
       columns(i,1:nrcol(i)) = (fscanf(fidmod,'%i',nrcol(i)))'; 
   %    disp(sprintf('Nrcol(i): ,length(col),sum(col): %i %i %i\n',nrcol(i),length(find(columns(i,:))),sum(columns(i,1:nrcol(i)))))
   end    
    cline = fgets(fidmod);     %#ok<NASGU> %SUN seems to need this extra line read?
    cline = fgets(fidmod);
   nexep  = sscanf(cline,'%*18c%i');  
   if nexep>0
      for i = 1:nexep
         l1(i)  = fscanf(fidmod,'%i',1);
         c1(i)  = fscanf(fidmod,'%i',1);
         l2(i)  = fscanf(fidmod,'%i',1);
         c2(i)  = fscanf(fidmod,'%i',1);
         pen(i) = fscanf(fidmod,'%f',1);
      end
   elseif nexep < 0     % DGM Nov 2006 - to support KKey's new exception method: blk1, blk2, penalty
        nexep = -nexep;
        for i = 1:nexep
            l1(i) = NaN;
            l2(i) = NaN;
            c1(i) = fscanf( fidmod, '%i', 1 );
            c2(i) = fscanf( fidmod, '%i', 1 );
            pen(i) = fscanf( fidmod, '%f', 1 );
        end
   else %assign l1 l2 c1 c2 pen = nan 
       l1  = nan;
       l2  = nan;
       c1  = nan;
       c2  = nan;
       pen = nan;
   end
   fclose(fidmod);
    return
    
%--------------------------------------------------------------------------
function[fidmsh,mshdesc,idx,nodey,nodez,nres,nfre,nexec,...
         fixres,freqs,meshy,meshz,nrc,nrl,mshres]...
          = readMESH(mshfile)
%--------------------------------------------------------------------------
% function readMESH.m
% Version 1.0, July 2001
% Written by: Kerry Key IGPP/SIO/UCSD
% 
% A matlab script to read an OCCAM2DMT MESH file
%
% Usage: 
%        [fidmsh,mshdesc,idx,nodey,nodez,nres,nfre,nexec,...
%         fixres,freqs,meshy,meshz,nrc,nrl,mshres]...
%          = readMESH(mshfile);
% 
% Inputs
%        mshfile: OCCAM2DMT MESH file name (string format)
%                 ***Can be FWD or INVERSE MESH
%
% Outputs: 
% 			fidmsh
% 			mshdesc
% 			idx
% 			nodey
% 			nodez
% 			nres
% 			nfre
% 			nexec
% 			fixres(nres)
% 			freqs(nfre)
% 			meshy(nodey-1)
% 			meshz(nodez-1)
% 			nrc
% 			nrl(nrc)
% 			mshres(4*(nodeyz-1),nodey-1)
%--------------------------------------------------------------------------

 
% OPEN AND READ IN MESH FILE:   

   fidmsh  = fopen(mshfile,'rt');
   mshdesc = fgets(fidmsh);           % mesh title/description
   idx     = fscanf(fidmsh,'%i',1);   % operation type, 0 to indicate Occam file
   nodey   = fscanf(fidmsh,'%i',1);   % num FE nodes in horizontal direction
   nodez   = fscanf(fidmsh,'%i',1);   % num FE nodes in vertical direction 
   nres    = fscanf(fidmsh,'%i',1);   % num fixed resistivities in mesh, up to 35
   nfre    = fscanf(fidmsh,'%i',1);   % num frequencies (ignored in inversion)
   nexec   = fscanf(fidmsh,'%i',1);   % nexec = 1 for fwd only, =2 for inversion, = 0 to print prameters and mesh
   
   if nres>0 && nres<=35
      fixres = fscanf(fidmsh,'%f',nres);  
   else
      fixres = nan ;
   end
   if nfre>0 
      freqs = fscanf(fidmsh,'%f',nfre); 
   else
      freqs = nan;
   end    
   
   meshy  = fscanf(fidmsh,'%f',nodey-1); % nodey-1 FE column widths, meters
   meshz  = fscanf(fidmsh,'%f',nodez-1); % nodez-1 FE layer heights, meters
   nrc    = fscanf(fidmsh,'%i',1);      % num stations, ignored by OCCAM
   if nrc>0
       nrl = fscanf(fidmsh,'%i',nrc);   %nrc number of nrl node numbers
   else 
       nrl = nan;
   end
   [mshres count] = fscanf(fidmsh,'%1s',[nodey-1 4*(nodez-1) ]);
   mshres = mshres';
   if count~=(nodey-1)*4*(nodez-1) 
       beep
       display('readMESH: Didn''t read correct number of mesh resistivity labels!!!')
   end
   
   
   fclose(fidmsh);

%--------------------------------------------------------------------------
function [YY, ZZ, restri] = makeELEMENTS(sumy,sumz,mshres)
%--------------------------------------------------------------------------    
% makeELEMENTS.m
%
% Written by:  Kerry Key
%              IGPP/SIO/UCSD
% Version:     1.0 August 2001
%
% About:   a function that takes mesh rox and column (y and z) distance vectors 
% along with an array of Wannamaker mesh resistivities and outputs the arrays
% that can plot up the mesh using the matlab intrinsic "patch.m"
% Returns RECTANGULAR ELEMENTS PLUS WITH SEAFLOOR BOUNDARY ELEMENTS AS TRIANGLES.
% Current removes rectanglur elements along seafloor, so no duplicates. 
%
% Inputs: sumy: cumulative sum of y-column increments (+ to right)
%         sumz: cumulative sum of z-column increments (- down)
%         mshres: array of mesh resistivities [4*(sumz-1) by sumy-1]
%                 -4 Rows per each row of sumz for triangles                
%--------------------------------------------------------------------------

   nodey = length(sumy);
   nodez = length(sumz);
% Need 4 coordinates for polygon vertices:
%  MAKE REPEATING IND.
   E1 = repmat([1 2 2 1]',1,nodey-1);
   E2 = repmat(1:nodey-1,4,1);
   Y  = E1+E2-1;
   Y  = repmat(Y,1,nodez-1);
   Z  = repmat([2 2 1 1]',1,(nodez-1))+repmat(0:nodez-2,4,1);
   Z  = repmat (Z,nodey-1,1);
   Z  = reshape(Z,4,(nodez-1)*(nodey-1));
   Z = fliplr(Z);
% NOW ASSIGN sum values:
   YY = sumy(Y);
   ZZ = sumz(Z);
% GET MESH ALPHANUMERICS INTO RES. TABLE INDICIES  (1,9,A-Z,? --> 1:35 + ?)
   temp = double(mshres);                % convert char array to double
   [I]  = find(temp<=57 & temp>= 48);    % mshres is 1-9 characters
      if ~isempty(I);
          temp(I)=temp(I)-48;
      end
   [I]  = find(temp>=65 & temp<=89);     % mshres is A-Y characters
      if ~isempty(I);
         temp(I) = temp(I)-55;
      end
   [Inan]  = find(temp==90);     %make the sea NAN
      if ~isempty(Inan);
         temp(Inan) = 999;
      end
   [I] = find(temp==63);                 % mshres is ?, free parameter
      if ~isempty(I);
          temp(I)= 63;
      end;       %not sure what to do with ? yet.
   lt  = size(temp);
   lt=lt(1);
% FIND TRAINGULAR ELEMENTS:
   [I,J]=find((temp(1:4:lt,:)~=temp(2:4:lt,:) | temp(1:4:lt,:)~=temp(3:4:lt,:) ...
              | temp(1:4:lt,:)~=temp(4:4:lt,:)) & temp(3:4:lt,:)~=999);
          
% reset 999 (sewater resistivities to nan (can't use find w/nan)
     temp(Inan)=nan;
%
     restri = reshape(temp(1:4:lt,:)',1,(nodey-1)*(nodez-1));
% APPEND ON TRIANGULAR ELEMENTS, ADDS 4 NEW COLUMN IND./TRI ELEMENT:
if ~isempty(I) && ~isempty(J)
   TRIEL  = sort(sub2ind(size(mshres'),J,I));
   % GET RID OF SEAFLOOR RECTANGLES, WILL USE TRIANGLES ALONG BOUNDARY
   INDX   = ones(1,length(YY));
   INDX(TRIEL) = 0;
   INDX = logical(INDX);
   mYY    = (YY(1,TRIEL)+YY(2,TRIEL))/2;
   YYT1   = [YY(1,TRIEL)' YY(2,TRIEL)' mYY' YY(1,TRIEL)']';
   YYT4   = [YY(2,TRIEL)' YY(3,TRIEL)' mYY' YY(2,TRIEL)']';
   YYT3   = [YY(3,TRIEL)' YY(4,TRIEL)' mYY' YY(3,TRIEL)']';
   YYT2   = [YY(4,TRIEL)' YY(1,TRIEL)' mYY' YY(4,TRIEL)']';
   YY     = [YY(:,INDX), YYT1, YYT2, YYT3, YYT4];

   mZZ    = (ZZ(1,TRIEL)+ZZ(3,TRIEL))/2;
   ZZT1   = [ZZ(1,TRIEL)' ZZ(2,TRIEL)' mZZ' ZZ(1,TRIEL)']';
   ZZT4   = [ZZ(2,TRIEL)' ZZ(3,TRIEL)' mZZ' ZZ(2,TRIEL)']';
   ZZT3   = [ZZ(3,TRIEL)' ZZ(4,TRIEL)' mZZ' ZZ(3,TRIEL)']';
   ZZT2   = [ZZ(4,TRIEL)' ZZ(1,TRIEL)' mZZ' ZZ(4,TRIEL)']';
   ZZ     = [ZZ(:,INDX), ZZT1, ZZT2, ZZT3, ZZT4];

   I      = (I-1)*4+1;
   TRIEL  = sort(sub2ind(size(mshres'),J,I));
   IND1   = TRIEL;
   IND2   = IND1+nodey-1;
   IND3   = IND1+2*(nodey-1);
   IND4   = IND1+3*(nodey-1);
   temp2  = temp';
 
   restri = [restri(INDX), temp2(IND1)', temp2(IND2)', ...
              temp2(IND3)',temp2(IND4)' ];
end
%--------------------------------------------------------------------------
function [YLy,YLz,ZLy,ZLz ] = meshlines(sumy,sumz)
%--------------------------------------------------------------------------
%    function meshlines.m
%
% Written by:  Kerry Key
%              IGPP/SIO/UCSD
% Version:     1.0 August 2001
%
% About:  makes lines for plotting PW2D MESH grid lines
%
% Usage: [YLy YLz ZLy ZLz ] = meshlines(sumy,sumz)
%
% Plotting the results: h_mesh = plot(YLy,YLz,'k-',ZLy,ZLz,'k-'); 
% Turn mesh off:        set(h_mesh,'visible','off');
% Turn on mesh:         set(h_mesh,'visible','on');
%--------------------------------------------------------------------------
   


% HORIZONTAL LINES:
   lsumy=length(sumy);
   lsumz=length(sumz);
   YLy=[sumy(1).*ones(lsumz,1) ,sumy(lsumy).*ones(lsumz,1)]';
   YLz=[sumz sumz]';
% VERTICAL LINES:
   ZLz=[sumz(1).*ones(lsumy,1) ,sumz(lsumz).*ones(lsumy,1)]';
   ZLy=[sumy sumy]';
 
 % h_mesh=plot(YLy,YLz,'k-',ZLy,ZLz,'k-'); axis ij
 % to turn mesh off: set(h_mesh,'visible','off');
 % to trun on mesh:  set(h_mesh,'visible','on');
 
%-------------------------------------------------------------------------- 
 function [REGY, REGZ, freebrick, seabrick, fixedbrick] = ...
           assembleREG(nlayer,nrcol,columns,irz,mshres,sumy,sumz)
%--------------------------------------------------------------------------
%              function assembleREG.m
%
% Written by:  Kerry Key
%              IGPP/SIO/UCSD
%
% Version:     1.0 August 2001
%
% About:       A function to assemble an OCCAM2DMT regularization grid,. 
%              including consolidated bricks lower down in grid.
%
% Usage:       [REGY, REGZ, freebrick, seabrick, fixedbrick] =
%              assembleREG(nlayer,nrcol,columns,irz,mshres,sumy,sumz);
% 
% Inputs:      nlayer, nrcol, columns, irz, mshres, sumy, sumz 
%
% Outputs:     REGY:       [4 x NBRICKS] array of BRICK Y positions
%              REGZ:       [4 x NBRICKS] array of BRICK Z positions
%              freebrick:  array freebrick indices in REGY and REGZ
%                          (bricks have free inversion parameters)  
%              fixedbrick: array fixedbrick indices in REGY and REGZ 
%                          (bricks are all fixed structure)
%              seabrick:   array seabrick indices in REGY and REGZ 
%                          (bricks are all seawater)
%--------------------------------------------------------------------------

% MAKE POSITION INDICES FOR GENERALIZED REG GRID 
%--------------------------------------------------------------------------
 % FORM REGY: ARRAY OF BRICK CORNER Y LOCATIONS   
   REGY = [];
   for i = 1:nlayer
       sumcols = [1; cumsum(columns(i,1:nrcol(i)))'+1]';
       regy = sumy(sumcols);
       E1 = repmat([1 2 2 1]',1,nrcol(i));
       E2 = repmat(1:nrcol(i),4,1); 
       Y  = E1+E2-1;
       REGY=[REGY regy(Y)];
   end
 % FORM REGZ: ARRAY OF BRICK CORNER Z LOCATIONS   
 % REG. ROW POSTION VECTOR  
   regz=sumz([1 cumsum(irz)+1]);    % distances in z, all columns the same
 % FORM BRICK DEPTH INDICES (1 PER LAYER): 
   Z  = repmat([1 1 2 2]',1,nlayer)+repmat(0:nlayer-1,4,1);
 % NOW REPLICATE nrcol(i) TIMES PER LAYER i
   Z2 = [];
   for i = 1:nlayer                        
      Z2 =[Z2 repmat(Z(:,i),1,nrcol(i))];  
   end
   REGZ=regz(Z2);  % ARRAY OF DEPTHS (m) [4 X #BRICKS]     
 %--------------------------------------------------------------------------
 % SEARCH FOR FREE BRICKS
 %--------------------------------------------------------------------------
 % SEARCH ORDER:  if freebrick, elseif seabrick, elseif fixedbrick

 %  Search bricks starting at top left corner (as res. values are given in 
 %  iteration files.
  seabrick   = [];
  fixedbrick = [];
   inz1    = 1;
   ibrknum = 0;
   freecnt = 1;
   fixcnt  = 1;
   seacnt  = 1;
   for i = 1:nlayer
       inz2 = inz1 + irz(i) - 1;
       iny1 = 1;
% for every brick in that layer 
       for j = 1:nrcol(i)
          ibrknum = ibrknum+1;
% find the right edge of brick
          iny2 = iny1 + abs(columns(i,j)) - 1;
% test to see if the brick is free 
          isfree = 0;
          issea  = 0;
          isfix  = 0;
          for ii = inz1:inz2
             for jj = iny1: iny2 
                if mshres(4*(ii-1)+3,jj) =='?'  %subtriangle #3 is bottom of element 
                   isfree = 1;
                   break;
                elseif mshres(4*(ii-1)+3,jj) =='Z' 
                   issea = 1;
                elseif mshres(4*(ii-1)+3,jj) ~='?' && ...
                       mshres(4*(ii-1)+3,jj) ~='Z'
                   isfix = 1;
                end   
             end
             if isfree == 1
                freebrick(freecnt) = ibrknum;
                freecnt = freecnt + 1;
                break
             end
          end
          if isfree~=1 && issea == 1
             seabrick(seacnt) = ibrknum;
             seacnt = seacnt + 1;
          elseif isfree~=1 && isfix == 1
             fixedbrick(fixcnt) = ibrknum;
             fixcnt = fixcnt + 1; 
          end
          iny1 = iny2 + 1;
       end   % j loop
       inz1 = inz2 + 1;
    end    % i loop
    if isempty(seabrick)
        seabrick = nan;
    end
    if isempty(fixedbrick)
        fixedbrick = nan;
    end    

%--------------------------------------------------------------------------
function [datfrm,datdes,nsites,sites,offsts,nfre,freqs, ...
        ndblks,DATA]= readDATA(datfile)
%--------------------------------------------------------------------------
% function readDATA.m
% Version 1.0, July 2001
% Written by: Kerry Key IGPP/SIO/UCSD
% 
% A matlab script to read an OCCAM2DMT DATA file
%
% Usage: 
%        [datfrm,datdes,nsites,sites,offsts,nfre,freqs, ...
%         ndblks,DATA]= readDATA(datfile);
%
% Inputs
%        datfile: OCCAM2DMT DATA file name (string format)
%
% Outputs: 
% 		datfrm
% 		datdes
% 		nsites
% 		sites(nsites)
% 		offsts(nsites)
% 		nfre:            Also can be found in MESH FILE, if running FWD PW2D 
% 		freqs(nfre):     Also can be found in  MESH FILE, if running FWD PW2D
% 		ndblks:          Number of data blocks
% 		DATA(ndblks,5):  Ndblks by [SITE  FREQ#  TYPE  DATUM  ERR]
%
%--------------------------------------------------------------------------


% OPEN AND READ IN DATA FILE (to get site offsets):
   fiddat = fopen(datfile,'r');
    cline  = fgets(fiddat);  
   datfrm = sscanf(cline,'FORMAT:           %80c');
   if  strcmp(deblank(datfrm),'OCCAM2DMTDATA_2.0')
      ndp = 4;
   else
       ndp =3;
   end
   
    cline  = fgets(fiddat); 
   datdes = sscanf(cline,'DESCRIPTION:      %80c');
    cline  = fgets(fiddat); 
   nsites = sscanf(cline,'SITES:            %f');
   for i  = 1:nsites
       cline    = fgets(fiddat);
      sites(i) = cellstr(sscanf(cline,'%s'));
   end
   string = fgets(fiddat);   %#ok<NASGU>
   offsts = fscanf(fiddat,'%f',nsites);
    cline  = fgets(fiddat);  %#ok<NASGU>
    cline  = fgetl(fiddat); 
   nfre   = sscanf(cline,'FREQUENCIES:%f');
   freqs  = fscanf(fiddat,'%f',nfre);
    cline  = fgets(fiddat); %#ok<NASGU>
    cline  = fgets(fiddat);
   ndblks = sscanf(cline,'DATA BLOCKS:%f');
    string = fgets(fiddat); %#ok<NASGU>
   DATA   = fscanf(fiddat,'%f %f %f %f %f',[ndp+2 ndblks]);
   DATA=DATA';
   fclose(fiddat);