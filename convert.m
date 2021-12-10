function convert(data_dir)

if ~exist('data_dir','var') || isempty(data_dir)
    data_dir='data';
end

years=get_subdir_names(data_dir);

for i=1:numel(years)
    models=get_subdir_names(fullfile(data_dir,years{i}));
    for j=1:numel(models)
        months=get_subdir_names(fullfile(data_dir,years{i},models{j}));
        for k=1:numel(months)
            files=get_subdir_names(fullfile(data_dir,years{i},models{j},months{k},'*.180'));
            for l=1:numel(files)
                file_now=fullfile(data_dir,years{i},models{j},months{k},files{l});
                if isempty(dir([file_now,'.mat']))
                    if isempty(dir(file_now))
                        error([mfilename,': cannot find model ',files{l},' in either ICGEM or mat formats.'])
                    else
                        disp(['Converting the following model to mat format: ',files{l}])
                        mod_save(mod_dyn_input(file_now,'struct'),file_now,'mat')
                    end
                else
                    disp(['The following model is already in mat format: ',files{l}])
                end
                if ~isempty(dir(file_now))
                    delete(file_now)
                end
            end
        end
    end
end

end

%% aux functions

function out=get_subdir_names(dir_name)
   d=dir(dir_name);
   out=remove_dot_dirs({d.name});
end

function out=remove_dot_dirs(dir_names)
    c=0;
    out=cell(0);
    for i=1:numel(dir_names(:))
        if dir_names{i}(1)=='.'
            %do nothing
        else
            c=c+1;
            out(c)=dir_names(i);
        end
    end
end

function out=mod_dyn_input(in,mode)
  % OUT=MOD_DYN_INPUT(IN,MODE) reads input IN and depending on it's type,
  % returns OUT
  %
  %   Decision tree:
  %
  %   type of IN:
  %
  %   - char: input is a filename (complete with paths) to supported files
  %   types (see below).
  %   - structure: input is a <mod> structures (see mod_struct2data).
  %   - 4-by-n arrays: input is a matrice with columns degree, order, C
  %   coefficient and S coefficient.
  %   - cell array: elements can be of any of the above types.
  %
  %   Output OUT is a cell array (with the same size as IN if the latter is
  %   a cell array) containing:
  %   - a matrix if IN is:
  %       - char and refers to a 'geoidDeg', 'geoidGrid', 'eqwDeg' or 'eqwGrid' file
  %       - 3-by-n matrix (lat,long,geoid)
  %       - 2-by-n matrix (degre,degree amplitudes)
  %   - a <mod> structure if IN is:
  %       - char and refers to a 'mod' file
  %       - a mod strucutre
  %       - a 4-by-n array (degree,order,C coefficient, S coefficient)
  %
  %   Supported files types:
  %       mod, geoidDeg,  geoidGrid, eqwDeg, eqwGrid
  %
  %   The optional input MODE allows this routine to check if a specific type
  %   of output is returned. Default is blank, meaning that no check is done.
  %   Supported are:
  %
  %   'scalar cell' : OUT is a cell array with one entry
  %   'non-cell'    : OUT is either a <mod> structure or a matrix (containing
  %                   geoidDeg,geoidGrid,eqwDeg or eqwGrid data)
  %   'struct'      : OUT is a <mod> structured (automatic conversion from
  %                   scalar cell)

  % Created by J.Encarnacao <J.deTeixeiradaEncarnacao@tudelft.nl>

  if ~exist('mode','var')
      mode='';
  end

  %determining input data type
  switch class(in)
      case 'char'
          if contains(in,'*') || ...
             contains(in,'|') || ...
             contains(in,'[') || ...
             contains(in,']') || ...
             contains(in,'^') || ...
             contains(in,'$')
              %wildcard filename
              out=mod_dyn_input(filename_wildcard(in));
          else
              if isempty(in)
                  error([mfilename,': cannot deal with empty filenames,'])
              end
              %getting file extension
              [~,~,e]=fileparts(in);
              %loading data from a file
              switch lower(e)
                  case {'.mod','.txt'}
                      out={mod_load(in)};
                      %checking for unknown GM and R
                      if out{1}.GM==-1 || out{1}.R==-1
                          %GM and R are the default ones
                          out=mod_data2struct(out);
                      end
                  case {'.gfc','.180'}
                      out={mod_load_ICGEM(in)};
                      %checking for unknown GM and R
                      if out{1}.GM==-1 || out{1}.R==-1
                          %GM and R are the default ones
                          out=mod_data2struct(out);
                      end
                  case {'.mat'}
                      load(in,'mod');
                      if ~iscell(mod)
                          out={mod};
                      end
                  % case {'.gsm'}
                  %     out={mod_load_GSM(in)};
                  %     %checking for unknown GM and R
                  %     if out{1}.GM==-1 || out{1}.R==-1
                  %         %GM and R are the default ones
                  %         out=mod_data2struct(out);
                  %     end
                  case {'.geoiddeg','.eqwdeg','.geoidgrid','.eqwgrid'}
                      out=file_load_datasets({in});
                  otherwise
                      error([mfilename,': unsupported filetype: ',e])
              end
          end
      case 'double'
          %determining what is the type of this data set
          switch size(in,2)
              case {2,3} %grid or deg file
                  out={in};
              case 4
                  %GM and R are the default ones
                  out={mod_data2struct(in)};
              otherwise
                  error([mfilename,': the size of the input data is ilegal: ',num2str(size(in,2))])
          end
      case 'cell'
          %going recursively through each elements of this cell array
           out=f_vector_call(@mod_dyn_input,1,in,'non-cell');
  %         out=cell(size(in));
  %         for i=1:length(in(:))
  %             out(i)=mod_dyn_input(in{i});
  %         end
      case 'struct'
          %input checking
          mod_struct2data(in);
          out={in};
      otherwise
          error([mfilename,': ilegal class for input IN: ',class(in),'.'])
  end

  %bug trap
  if ~iscell(out)
      error([mfilename,':BUG: expecting this routine to return always cells; now returning ',class(out)])
  end

  %delivering requested types
  switch lower(mode)
      case('scalar cell')
          if max(size(out)) > 1
              error([mfilename,': input IN must have at most 1 entry.'])
          end
          if ~iscell(out)
              out={out};
          end
      case('struct')
          if max(size(out)) > 1
              error([mfilename,': input IN must have at most 1 entry.'])
          end
          %checking if this is the correct structure
          mod_struct2data(out{1});
          %propagating
          out=out{1};
      case ('non-cell')
          if max(size(out)) > 1
              error([mfilename,': input IN must have at most 1 entry.'])
          end
          %propagating
          out=out{1};
      case('')
          %do nothing
      otherwise
          error([mfilename,': unknown <mode> ',mode])
  end
end

function mod=mod_load(filenames)

  if ~iscell(filenames) && ischar(filenames)
      filenames={filenames};
  elseif ~iscellstr(filenames)
      error([mfilename,': input <filenames> must be a string or cell of strings.'])
  end

  %loading data (always reading the original file, there's no need for saving it in mat format)
  models=file_load_datasets(filenames);

  %loading constants
  GM=cell(1,length(filenames(:)));
   R=cell(1,length(filenames(:)));
  GM(:)={-1};
   R(:)={-1};
  %going through all files
  for i=1:length(filenames(:))
      %opening this file
      fid = fopen_disp(filenames{i});
      %searching for ID strings
      while GM{i}==-1 || R{i}==-1
          c=textscan(fid,'%s',1);
          %checking for end-of-file
          if isempty(c{1}) || feof(fid)
              disp(['WARNING: Could not find GM (',num2str(GM{i}),') and/or  R (',num2str(R{i}),').'])
              break
          end
          %checking if this is one of the ID strings
          switch c{1}{1}
              case 'GM:'
                  c=textscan(fid,'%f',1);
                  GM{i}=c{1}(1);
              case 'R:'
                  c=textscan(fid,'%f',1);
                  R{i}=c{1}(1);
          end
      end
      %closing file
      fclose(fid);
  end

  %checking all models
  mod=mod_data2struct(models,GM,R);
end

function [mod,mod_error]=mod_load_ICGEM(filenames)

  if ~iscell(filenames) && ischar(filenames)
      filenames={filenames};
  elseif ~iscellstr(filenames)
      error([mfilename,': input <filenames> must be a string or cell of strings.'])
  end

  mod=cell(size(filenames));
  mod_error=cell(size(filenames));

  %going through all files
  for i=1:length(filenames(:))

      %loading the data
      [cnm,snm,ecnm,esnm,header]=load_icgem(filenames{i});

      %converting to internal format
      mod{i}=mod_matrix2list(cnm,snm,header.earth_gravity_constant,header.radius);

      if ~isempty(ecnm) && ~isempty(esnm)
          mod_error{i}=mod_matrix2list(ecnm,esnm,header.earth_gravity_constant,header.radius);
      end

  end

  if length(mod(:)) == 1
      mod=mod{1};
      if length(mod_error(:)) == 1
          mod_error=mod_error{1};
      end
  end
end

%% auxiliar function

function [cnm,snm,ecnm,esnm,header,modelname,n_t0,n_trnd,n_acos,n_asin,cnm_t0,cnm_trnd,snm_trnd,cnm_acos,snm_acos,cnm_asin,snm_asin]=load_icgem(filename)

%This function is an adaptation of icgem2mat.m from rotating_3d_globe, by
%Ales Bezdek, which can be found at:
%
%http://www.asu.cas.cz/~bezdek/vyzkum/rotating_3d_globe/
%
%The original header is transcribed below.
%
%J.Encarnacao (j.g.deteixeiradaencarnacao@tudelft.nl) 11/2013

% ICGEM2MAT   Reads geopotential coefficients from an ICGEM file and saves them in a mat file.
%
% Usage:
%
%       icgem2mat
%
% finds all the ICGEM files (*.gfc) in the current directory,
% reads the geopotential coefficients, transforms them into Matlab variables:
%       header...structure with Icgem header information
%       cnm(n+1,m+1), snm(n+1,m+1)...harmonic coefficients C(n,m), S(n,m)
%
% The new mat file with the same name is moved into 'data_icgem' subdirectory;
% the original gfc file is moved into 'data_icgem/gfc/' subdirectory.
%
% Add the 'data_icgem' folder into your Matlab path.
% The model coefficients are then loaded by typing, e.g.:
%
%          load egm2008
%
% To display the C(2,0) zonal term type
%
%          cnm(3,1)
%
%
% See also compute_geopot_grids

% Ales Bezdek, bezdek@asu.cas.cz, 11/2012

% clear
% NMAX=360;
% NMAX=1e100;  %it is possible to limit the maximum degree read from the gfc file
% adr_data='./';
% adr_kam='./data_icgem/';

% seznam_soub=dir(adr_data);
% soub={seznam_soub.name};   %cell with filenames
% for i=1:length(soub)
%    jm=soub{i};
%    if length(jm)>4 && strcmpi(jm(end-3:end),'.gfc')
%       soub1=jm(1:end-4);
%       fprintf('Gfc file processed: %s\n',soub1);
%       filename=[adr_data soub1 '.gfc'];

      % Read header
      fid=fopen(filename);
      modelname=''; GM=0; ae=0; Lmax=0; errors=''; norm=''; tide='';
      product_type='';model_content='';

      s=fgets(fid);
      while(strncmp(s, 'end_of_head', 11) == 0 && sum(s)>=0)
         if (keyword_search(s, 'product_type'          )), product_type =strtrim(s(13:end)); end;
         if (keyword_search(s, 'modelname'             )), modelname    =strtrim(s(10:end)); end;
         if (keyword_search(s, 'model_content'         )), model_content=strtrim(s(14:end)); end;
         if (keyword_search(s, 'errors'                )), errors       =strtrim(s(7:end)); end;
         if (keyword_search(s, 'norm'                  )), norm         =strtrim(s(5:end)); end;
         if (keyword_search(s, 'tide_system'           )), tide         =strtrim(s(12:end)); end;
         if (keyword_search(s, 'earth_gravity_constant')), GM           =str2double(strrep(s(23:end),'D','e')); end;
         if (keyword_search(s, 'radius'                )), ae           =str2double(s(7:end)); end;
         if (keyword_search(s, 'max_degree'            )), Lmax         =str2double(s(11:end)); end;
         s=fgets(fid);
      end
      if sum(s)<0
         error_ab('Problem with reading the gfc file.')
      end

      %J.Encarnacao: added the filename
      header=struct(...
          'product_type',           product_type,...
          'modelname',              modelname,...
          'model_content',          model_content,...
          'earth_gravity_constant', GM,...
          'radius',                 ae,...
          'max_degree',             Lmax,...
          'errors',                 errors,...
          'norm',                   norm,...
          'tide_system',            tide,...
          'filename',               filename...
      );

      % read coefficients
      cnm=zeros(Lmax+1);
      snm=zeros(Lmax+1);
      ecnm=zeros(Lmax+1);
      esnm=zeros(Lmax+1);

      i_t0=0;
      i_trnd=0; %pocet clenu s trendem
      i_acos=0; %pocet clenu
      i_asin=0; %pocet clenu
      i_gfc=0;

     cnm_t0=[]; cnm_trnd=[]; snm_trnd=[]; cnm_acos=[]; snm_acos=[]; cnm_asin=[]; snm_asin=[];

      s=fgets(fid);
      while (s>=0)

        %skip empty lines
        if numel(s)<5
          s=fgets(fid);
          continue
        end

         x=str2num(s(5:end));
         n=x(1)+1;
         m=x(2)+1;
%          if n>NMAX || m>NMAX
%             s=fgets(fid);
%             continue;
%          end
%         disp(s(1:4))
         if strcmp(s(1:4),'gfct')
            if isempty(cnm_t0)
               [~, result] =system(['grep -c gfct ' filename]);
               i1=str2double(result);
               if i1==0; error_ab('Problem with t0'); end
               cnm_t0=zeros(i1,3);
            end
            i_t0=i_t0+1;
            cnm(n,m)=x(3);
            snm(n,m)=x(4);
%            if strcmp(header.errors, 'formal') || strcmp(header.errors,'calibrated')
               [yr,mn,dy]=ymd2cal(x(end)/1e4);
               yrd=jd2yr(cal2jd(yr,mn,dy));
               cnm_t0(i_t0,:)=[n m yrd];
%             elseif strcmp(header.errors,'calibrated_and_formal')
%             elseif strcmp(header.errors,'no')
%             end
            if (strcmp(header.errors, 'formal') || strcmp(header.errors,'calibrated') || strcmp(header.errors,'calibrated_and_formal')),
               ecnm(n,m)=x(5);
               esnm(n,m)=x(6);
            end
         elseif strcmp(s(1:3),'gfc')
            cnm(n,m)=x(3);
            snm(n,m)=x(4);
            if (strcmp(header.errors, 'formal') || strcmp(header.errors,'calibrated') || strcmp(header.errors,'calibrated_and_formal')),
               ecnm(n,m)=x(5);
               esnm(n,m)=x(6);
            end
            i_gfc=i_gfc+1;
         elseif strcmp(s(1:4),'trnd') || strcmp(s(1:3),'dot')
            if isempty(cnm_trnd)
               [~, result] =system(['grep -c trnd ' filename]);
               i1=str2double(result);
               if i1==0; [~, result] =system(['grep -c dot ' filename]); i1=str2num(result); end
               if i1==0; error_ab('Problem with trnd'); end
               cnm_trnd=zeros(i1,3); snm_trnd=cnm_trnd;
            end
            i_trnd=i_trnd+1;
            cnm_trnd(i_trnd,:)=[n m x(3)];
            snm_trnd(i_trnd,:)=[n m x(4)];
         elseif strcmp(s(1:4),'acos')
            if isempty(cnm_acos)
               [~, result] =system(['grep -c acos ' filename]);
               i1=str2double(result);
               if i1==0; error_ab('Problem with acos'); end
               cnm_acos=zeros(i1,4); snm_acos=cnm_acos;
            end
            i_acos=i_acos+1;
            cnm_acos(i_acos,:)=[n m x(3) x(end)];
            snm_acos(i_acos,:)=[n m x(4) x(end)];
         elseif strcmp(s(1:4),'asin')
            if isempty(cnm_asin)
               [~, result] =system(['grep -c asin ' filename]);
               i1=str2double(result);
               if i1==0; error_ab('Problem with asin'); end
               cnm_asin=zeros(i1,4); snm_asin=cnm_asin;
            end
            i_asin=i_asin+1;
            cnm_asin(i_asin,:)=[n m x(3) x(end)];
            snm_asin(i_asin,:)=[n m x(4) x(end)];
         else
            error_ab('A problem occured in gfc data.');
         end
         s=fgets(fid);
      end
      fclose(fid);

      %J.Encarnacao: handle the tide system
      switch header.tide_system
      case 'zero_tide'
        %do nothing, this is the default
      case {'free_tide','tide_free'}
        cnm(3,1)=cnm(3,1)-4.173e-9;
        header.tide_system='zero_tide';
      case 'mean_tide'
        %http://mitgcm.org/~mlosch/geoidcookbook/node9.html
        cnm(3,1)=cnm(3,1)+1.39e-8;
        header.tide_system='zero_tide';
      otherwise
        %the tide system is not documented, so make some assumptions
        switch header.modelname
        case 'GROOPS'
          %do nothing, Norber Zehentner's solutions are zero tide
        case 'Improved Mass Transport Model'
          %do nothing, i couldn't find any information, so I assume it's zero tide
        otherwise
          error([mfilename,': unknown tide system ''',header.tide_system,'''.'])
        end
      end

%       modelname=header.modelname;

      %it is possible to limit the maximum degree read from the gfc file
      n_gfc=i_gfc; n_t0=i_t0; n_trnd=i_trnd; n_acos=i_acos; n_asin=i_asin;
      if n_t0; cnm_t0=cnm_t0(1:n_t0,:); end
      if n_trnd; cnm_trnd=cnm_trnd(1:n_trnd,:); snm_trnd=snm_trnd(1:n_trnd,:); end
      if n_acos; cnm_acos=cnm_acos(1:n_acos,:); snm_acos=snm_acos(1:n_acos,:); end
      if n_asin; cnm_asin=cnm_asin(1:n_asin,:); snm_asin=snm_asin(1:n_asin,:); end
      if n_t0~=n_trnd || n_acos~=n_asin
         error_ab('Problem with numbers of TVG terms.');
      end
%       fprintf('   gfc terms: %d, gfct: %d, trnd: %d, acos: %d, asin: %d\n',[n_gfc n_t0 n_trnd n_acos n_asin]);

%       if ~exist(adr_kam,'file'), mkdir(adr_kam); end
%       if ~exist([adr_kam 'gfc'],'file'), mkdir([adr_kam 'gfc']); end
%       eval(sprintf('save %s%s.mat cnm snm ecnm esnm header modelname n_t0 n_trnd n_acos n_asin cnm_t0 cnm_trnd snm_trnd cnm_acos snm_acos cnm_asin snm_asin;',adr_kam,soub1));
%       movefile([adr_data soub1 '.gfc'],[adr_kam 'gfc']);
%       fprintf('  Resulting file %s.mat was moved into folder: %s\n',soub1,adr_kam);
%       fprintf('  Original file  %s.gfc was moved into folder: %sgfc\n',soub1,adr_kam);
%    end
end


function out=keyword_search(line,keyword)
    out=strncmp(line,       keyword,         length(keyword)) || ...
        strncmp(line,strrep(keyword,' ','_'),length(keyword));
end

function mod_save(mod,filename,mode)
  % MOD_SAVE(MOD,FILENAME) saves MOD structure in FILE.

  % Created by J.Encarnacao <J.G.deTeixeiradaEncarnacao@tudelft.nl>

  if ~exist('mode','var') || isempty(mode)
      mode='mod';
  end

  if all([iscell(mod) iscell(filename)])
      if any(size(mod) ~= size(filename))
         error([mfilename,': if <mod> and <filename> are cells, they need to have the same size'])
      end

      f_vector_call(@mod_save,[1,2],mod,filename,mode);
      return
  elseif xor(iscell(mod),iscell(filename))
      error([mfilename,': if either <mod> or <filename> is cell array, the other must also be cell array.'])
  end

  %need structure, not cell array
  mod=mod_dyn_input(mod,'struct');

  if (mod.mod(end,1)>9999)
      error([mfilename,': this routine only supports models up to degree and order 9999.'])
  end

  switch mode
      case 'mod'
          %opening output file
          fid=fopen(filename,'wt');

          %writting header
          fprintf(fid,'%s\n','#MOD-5.1 (Model of the Earth gravity field)');
          fprintf(fid,' %i %i | Lmin & Lmax\n',round(min(mod.mod(:,1))),round(max(mod.mod(:,1))));
          fprintf(fid,' %s\n',['File creator: matlab routine mod_save.m at ',datestr(now)]);
          fprintf(fid,'\n');
          fprintf(fid,'\n');
          fprintf(fid,'\n');
          fprintf(fid,'\n');
          fprintf(fid,'\n');
          fprintf(fid,'\n');
          fprintf(fid,'\n');
          fprintf(fid,'\n');
          fprintf(fid,'\n');
          fprintf(fid,'\n');
          fprintf(fid,'\n');
          fprintf(fid,'\n');
          fprintf(fid,'\n');
          fprintf(fid,'\n');
          fprintf(fid,' GM: %18.1f       m^3 s^-2; R:   %12.4f\n',mod.GM,mod.R);
          fprintf(fid,'\n');
          fprintf(fid,'   l    m             C                        S\n');

          %writting data
          for i=1:size(mod.mod,1)
              fprintf(fid,'%4d %4d % .16e % .16e\n',mod.mod(i,:));
          end

          %closing file
          fclose(fid);
      case 'mat'
          if ~strcmp(filename(end-3:end),'.mat')
              filename=[filename,'.mat'];
          end
          save(filename,'mod','-mat')
      otherwise
          error([mfilename,': unknown mode ',mode,'.'])
  end
end

function [o0,o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13,o14,o15,o16,o17,o18,o19]=f_vector_call(f_handle,cell_idx,varargin) %#ok<STOUT>
  % [O0,01,...O19] = F_VECTOR_CALL(FUN_HANLDE,CELL_IDX,VARARGIN) calls
  % FUNC_HANDLE several times, storing the outputs in outputs O0,O1, ..., O19.
  % The number of times that it is called depends on the size of the
  % CELL_IDX-th input (which must be a cell array) in varargin. In order for
  % this to work, VARARGIN must be compatible with the argument list of the
  % function specified in FUN_HANDLE
  %
  %   Example:
  %   % to plot several curves
  %   x=linspace(0,pi)
  %   y={sin(x),cos(x),sin(2*x)};
  %   fmt_str={'b','r','g'}
  %   %... using this function:
  %   f_vector_call(@plot,[2,3],x,y,fmt_str)
  %   %... is the same as:
  %   for i=1:length(y)
  %       plot(x,y{i},fmt_str{i})
  %       hold on %<f_vector_call> automatically calls <hold on>
  %   end
  %
  %   If there are outputs:
  %
  %   z=f_vector_call(@mean,1,y)
  %   %... is the same as:
  %   for i=1:length(y)
  %       z{i}=mean(y{i})
  %   end
  %
  %   A maximum of 19 outputs is supported.

  % Created by J.Encarnacao <J.G.deTeixeiradaEncarnacao@tudelft.nl>

  if ~exist(func2str(f_handle),'file') && ~exist(func2str(f_handle),'builtin')
      error([mfilename,': could not find function <',func2str(f_handle),'>.'])
  end
  if min(size(cell_idx)) ~= 1
      error([mfilename,': input <cell_idx> must be a vector of integer indexes.'])
  end

  n=size(varargin{cell_idx(1)});
  if min(n) == 0
     error([mfilename,': cell arguments cannot be empty.'])
  end
  if min(n) ~= 1
     error([mfilename,': cell arguments must be in cell vectors, not matrices (or higher dimension).'])
  end
  for i=1:length(cell_idx)
      if ~iscell(varargin{cell_idx(i)})
          error([mfilename,': input argument ',num2str(i),' of function ',func2str(f_handle),...
              ' must be a cell vector.'])
      end
      if any(size(varargin{cell_idx(i)}) ~= n)
          error([mfilename,': input argument ',num2str(i),' of function ',func2str(f_handle),...
              ' must have the same size as other inputs.'])
      end
  end

  %checking if too many output arguments are being asked
  if nargout > nargout(f_handle)
      error([mfilename,': function <',func2str(f_handle),'> has a max of ',num2str(nargout(f_handle)),...
          ' output arguments, while ',num2str(nargout),' are being requested.'])
  end

  argout=cell(max(n),nargout);

  %building argument list
  for i=1:max(n)
      argin=cell(1,length(varargin));
      for j=1:length(varargin)
          if any(j==cell_idx)
             argin(j)= varargin{j}(i);
          else
             argin(j)= varargin(j);
          end
      end
      if nargout == 0
          f_handle(argin{:});
      else
          [argout{i,:}]=f_handle(argin{:});
      end
      %in case this is a plotting routine
      if ~isempty(get(0,'CurrentFigure'))
          hold on
      end
  end

  %propagating outputs
  for i=1:nargout
      if n(1) == 1
          eval(['o',num2str(i-1),'=argout(:,i)'';'])
      else
          eval(['o',num2str(i-1),'=argout(:,i);'])
      end
  end
end

function out=mod_matrix2list(c,s,GM,R)

    %builds the model <out> out of cosine and sine coefficient matrix <c> and
    %<s>.

    if ~exist('GM','var') || isempty(GM)
        GM=[];
    end
    if ~exist('R','var') || isempty(R)
        R=[];
    end
    if size(c,1) ~= size(c,2)
        error([mfilename,': input matrices are not square.'])
    end
    if size(s,1) ~= size(c,1)
        error([mfilename,': input sine and consine matrices are not the same.'])
    end

    %filtering upper diag matrix
    c(triu_idx(c,1))=NaN;
    s(triu_idx(s,1))=NaN;

    %put it back to a list
    [ci,cj,cv]=matrix2idx(c);
    [si,sj,sv]=matrix2idx(s);

    %need to do some sorting
    [ci,idx]=sort(ci);cj=cj(idx);cv=cv(idx);
    [si,idx]=sort(si);sj=sj(idx);sv=sv(idx);

    %bug trap
    if any(ci ~= si) || any(cj ~= sj)
        error([mfilename,': bug trap: sorting of indeces produced indexes in different positions.'])
    end

    %outputs
    out=mod_data2struct([ci-1,cj-1,cv,sv],GM,R);

    %need to remove the sine coefficients of zero orders
    out.mod(out.mod(:,2)==0,4)=0;
end

function idx=triu_idx(X,k)

    %returns the indexes of the upper triangular part of <X>. <X> is unchanged,
    %it is only used to estabilsh the required matrix size.

    if ~exist('k','var') || isempty(k)
        k=0;
    end

    idx=1:numel(X);
    idxm=zeros(size(X));
    idxm(:)=idx;

    idx=(idxm==triu(idxm,k));
end

function [i,j,v] = matrix2idx(m)

    %NaN entries in matrix <m> are not part of [i,j,v]

    [jm,im]=meshgrid(1:size(m,1),1:size(m,2));

    i=im(:);
    j=jm(:);
    v=m(:);

    %filtering out NaNs
    idx=isnan(v);
    i(idx)=[];
    j(idx)=[];
    v(idx)=[];

end

function mod_struct=mod_data2struct(mod,GM,R)

%Joins the structure <mod_struct> from it's components <mod>, <GM> and
%<R>. If <mod> is already a structure, then it checks <mod%mod> for correct
%dimensions and sets <mod%GM> and <mod%R> equalt to <GM> and <R>, resp. It
%is then returned as <mod_struct>.
%
%<mod_struct.mod> is a spherical harmonic expansions of a gravity field.
%Column 1 has degree, column 2 has order, column 3 has the cosine terms and
%column 4 has sine terms.
%
%Use this routine for sanity checks for <mod>, <GM> and <R>.

    switch class(mod)
        case 'struct'
            %if <mod> is already a structure, then do some sanity checks
            mod_struct2data(mod);
            %if <GM> and <R> are give, then propagate them to the structure
            if exist('GM','var') && ~isempty(GM) && exist('R','var') && ~isempty(R)
                mod_struct=mod_data2struct(mod.mod,GM,R);
            else
                %checking for unknown GM and R
                if isempty(mod.GM) || mod.GM<0
                    mod.GM=mod_default_const;
                end
                if isempty(mod.R) || mod.R<0
                    [dummy,mod.R]=mod_default_const;
                end
                %propagate the whole structure
                mod_struct=mod;
            end
            %done here
            return

        case 'cell'
            if ~exist('GM','var') || isempty(GM) || ~exist('R','var') || isempty(R)
                %vector call
                mod_struct=function_vector_call(@mod_data2struct,1,mod);
            else
                %if mod is a cell array, then recursively go through all entries
                if length(mod(:)) == 1
                    %checking for scalar GM and R
                    if ~iscell(GM),tmp{1}=GM;GM=tmp; end
                    if ~iscell(R), tmp{1}=R;  R=tmp; end
                    %reducing scalar cells
                    mod_struct=mod_data2struct(mod{1},GM{1},R{1});
                else
                    %vector call
                    mod_struct=function_vector_call(@mod_data2struct,1:3,mod,GM,R);
                end
            end
            %done here
            return

        case 'double'

            if ~exist('GM','var') || isempty(GM)
                GM=mod_default_const;
            end
            if ~exist('R','var') || isempty(R)
                [dummy,R]=mod_default_const;
            end
            if ~isnumeric(mod) || ~isnumeric(GM) || ~isnumeric(R)
                error([mfilename,': inputs <mod>, <GM> and <R> must all be numeric.'])
            end
            if size(mod,2) ~= 4
                error([mfilename,': input <mod> must be a 4-by-n matrix.'])
            end

            mod_struct=struct('mod',mod,'GM',GM,'R',R);

        otherwise
            error([mfilename,': ilegal class for input <in>: ',class(in),'.'])
    end
end

function [mod,GM,R]=mod_struct2data(mod_struct)

%Splits the structure <mod_struct> into it's components <mod>, <GM> and
%<R>. If <mod_struct> is not a structure, then it checks for correct
%dimensions and returns it along with default values for <GM> and <R>.

    if iscell(mod_struct)

        if length(mod_struct(:)) == 1
            %reducing scalar cells, avoids returning 1-entry cell arrays
            [mod,GM,R]=mod_struct2data(mod_struct{1});
        else
            %vector call
            [mod,GM,R]=function_vector_call(@mod_struct2data,1,mod_struct);
        end

    %     mod=cell(size(mod_struct));
    %     GM=zeros(size(mod_struct));
    %     R=zeros(size(mod_struct));
    %     for i=1:length(mod_struct(:))
    %         [mod{i},GM(i),R(i)]=mod_struct2data(mod_struct{i});
    %     end

        return
    end

    if any(~isfield(mod_struct,{'mod','GM','R'}))
        error([mfilename,': this does not appear to be a <mod> structure.'])
    end
    if size(mod_struct.mod,2) ~= 4
        error([mfilename,': input <mod_struct%mod> must be a 4-by-n matrix.'])
    end

    mod=mod_struct.mod;
    GM=mod_struct.GM;
    R=mod_struct.R;

end