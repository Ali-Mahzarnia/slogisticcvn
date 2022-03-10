clear;
filename = "/Users/ali/Desktop/mar/mice/MasterSheet_Experiments2021.xlsx" ;
connectomes_folder = "/Users/ali/Desktop/mar/connectivitymatrixwithshortnames/";

connectomes_type = 'mat';
allfiles_data = readtable(filename,'Sheet','CVN_20abb15');



data_source = 2;

if data_source ==2
    outpath='/Users/ali/Desktop/mar/mice/cvn/results'
    if ~exist(outpath, 'dir')
       mkdir(outpath)
    end
end

%genotypekeySet = {'APOE33', 'APOE23', 'APOE34', 'APOE44'};
%genotypevalueSet = [3 3 4 4];
%geno = containers.Map(genotypekeySet,genotypevalueSet);

%sexkeySet = {'M', 'F'};
%sexvalueSet = [1 2];
%sex = containers.Map(sexkeySet,sexvalueSet);

treatment = {'treadmill', 'wheel_only', 'sedentary'};
%treatmentvals = [1 1 3]; %only seden and nonseden
treatmentvals = [1 2 0]; %all three groups
treat = containers.Map(treatment,treatmentvals);


response_table = table('Size',[size(allfiles_data,1) 2],'VariableNames',{'Treatment', 'DWI'}, 'VariableTypes',{'string','string'});
response_array_init = zeros([size(allfiles_data,1) 2]);
for i = 1:size(allfiles_data,1)
    if size(cell2mat(table2array(allfiles_data(i,"DWI"))),1)>0

        temp=table2array(allfiles_data(i,'DWI'));
        temp=regexp(temp,'\d+(\.)?(\d+)?','match');
    response_array_init(i,1) =str2double([temp{:}]);
    response_array_init(i,2) = treat(cell2mat(table2array(allfiles_data(i,'Treatment'))));
    response_table(i,1) = num2cell(response_array_init(i,1));
    response_table(i,2) = num2cell(response_array_init(i,2));
    end
end
response_table=rmmissing(response_table);


connectivity = zeros(332,332,size(response_table,1));
if strcmp(connectomes_type,'xlsx')
    getpath = join([connectomes_folder,'*connec*.xlsx'],"");
elseif strcmp(connectomes_type,'mat')
    getpath = join([connectomes_folder,'*.mat'],"");
end

files = dir(getpath);
subjlist = zeros(size(response_table,1),1);
notfoundlist = zeros(size(response_table,1),1);
j=1;
l=1;
for i = 1:size(response_table,1)
    subjname = response_table{i,1};
  % if  ~(subjname==3394)
    found = 0;
    for file = files'
        subj = strsplit(file.name,'_');
        subj = subj{1};
        subj = subj(2:(end-4));
        if  str2double(subjname) == str2double(subj)
            if strcmp(connectomes_type,'xlsx')
                csv = readtable(join([connectomes_folder,file.name],""));
                csv.Var1{84} = 'ctx-rh-insula';
                csv = removevars(csv,{'Var1'});
                csv = table2array(csv);
                connectivity(:,:,j) = csv;
                found = 1
                break
            elseif strcmp(connectomes_type,'mat')
                A = load(join([connectomes_folder,file.name],""));  
                connectivity(:,:,j) = A.connectivity;
                found = 1;
                break
            end
        end
    end
    if found==1
        %display('found '+ string(subjname))
        subjlist(j) = subjname;
        j = j + 1;
    else
        %display('did not find '+ string(subjname))
        notfoundlist(l) = subjname;
        l = l + 1;
    end
  % end 
end

subjlist = subjlist(subjlist ~= 0);
notfoundlist = notfoundlist(notfoundlist ~= 0);
connectivity = connectivity(:,:,1:size(subjlist,1));
response_array = zeros(size(subjlist,1),size(response_array_init,2));
i=1;
for k = 1:numel(response_array_init)
    if ismember(response_array_init(k),subjlist)
        response_array(i,:) = response_array_init(k,:);
        i=i+1;
    end
end

%subselect = '_genotype_4';
subselect = '';
%subselect = '_NCgt40';
%subselect = 'norisk';

if size(subselect,2)>0
    if contains(subselect,'genotype')
        if contains(subselect,'3')
            APOE3 = response_array(:,2)==3;
            connectivity = connectivity(:,:,APOE3);
            response_array = response_array(APOE3,:);
            subjlist = subjlist(APOE3);
        elseif contains(subselect,'4')
            APOE4 = response_array(:,2)==4;
            connectivity = connectivity(:,:,APOE4);
            response_array = response_array(APOE4,:);
            subjlist = subjlist(APOE4);
        end
    elseif contains(subselect,'_NCgt40')
        %NCgt40 = find((response_array(:,3)>40).*(response_array(:,5)<2)); % response(response(:,5)==2, 1); # better be nonzero or code crashes
        NCgt40 = find((response_array(:,3)>40)); % response(response(:,5)==2, 1); # better be nonzero or code crashes
        connectivity = connectivity(:,:,NCgt40);
        response_array =response_array(NCgt40,:);
        subjlist = subjlist(NCgt40);
    end
     elseif contains(subselect,'norisk')
        %NCgt40 = find((response_array(:,3)>40).*(response_array(:,5)<2)); % response(response(:,5)==2, 1); # better be nonzero or code crashes
        norisk = find((response_array(:,5)<2)); % response(response(:,5)==2, 1); # better be nonzero or code crashes
        connectivity = connectivity(:,:,norisk);
        response_array =response_array(norisk,:);
        subjlist = subjlist(norisk);
end




%{
% % % % two step regression to take out the age's effect

%{

%%%% standadrizing through connectivity
for i=2:size(connectivity,1)
    for j=1:(i-1) % only the lower traingle  without diagonal elements
            
        if std(connectivity(i,j,:))>0 
            tempmean=mean(connectivity(i,j,:));
            tempsd=std(connectivity(i,j,:));
           connectivity(i,j,:)=(connectivity(i,j,:)-tempmean)/tempsd;
            connectivity(j,i,:)= connectivity(i,j,:);
        end
    end   
end
%%%%
%}


% first take the effect of other risks from age: age~ sex + gene + risk factor ....

age=response_array(:,3); risks=response_array(:,[2 4 5]); risks= [ ones(size(risks, 1),1) risks ]; % intercept needs one in design matrix
[~,~,r]=regress(age,risks); 
age=r; % replace age by the residual of this regression

% second step


%{
A=randn(5,5)
At = A.';
m = tril(true(size(At)),-1);
v = At(m).'
BB=A*0
BB(m)=v
BB=BB.'
BB-A
%}


origtraingle=NaN( size(connectivity,3)   ,  (size(connectivity,1)*(size(connectivity,2)-1))/2    );

for i=1:size(connectivity,3)  %extract the lower triangle
        A=connectivity(:,:,i);
        At = A.';
        m = tril(true(size(At)),-1);
        v = At(m).';
        if std(v)==0; display(i); end % no column with all zero values
        origtraingle(i,:)=v;
end


for j=1:size(origtraingle,2)  % each connectivity is regressed to age in a ismple lin regression
[~,~,r]=regress(origtraingle(:,j),age);
origtraingle(:,j)=r;
end   

% put them back in
for i=1:size(connectivity,3) 
        A=connectivity(:,:,i);
        BB=A*0;
        BB(m)=origtraingle(i,:);
        BB=BB.';
        connectivity(:,:,i)=BB+BB.';
end

%%%%

%}
%}



origtraingle=NaN( size(connectivity,3)   ,  (size(connectivity,1)*(size(connectivity,2)-1))/2    );

for i=1:size(connectivity,3)  %extract the lower triangle
        A=connectivity(:,:,i);
        At = A.';
        m = tril(true(size(At)),-1);
        v = At(m).';
        if std(v)==0; display(i); end % no column with all zero values
        origtraingle(i,:)=v;
end

response_array_path = join([outpath 'response_array' subselect '.mat'],"");
response_table_path= join([outpath 'response_table' subselect '.mat'],""); 
connectomes_path = join([outpath 'connectivity_all' '_ADDecode' '_Dipy' subselect '.mat'],"");
trainglepath = join([outpath 'traingle' subselect '.mat'],"");


save(response_array_path, 'response_array');
save(response_table_path, 'response_table');
save(connectomes_path, 'connectivity', 'subjlist');
save(trainglepath, 'origtraingle')


%{
%hstogram of ages with edge of bins:
edge1=35;
edge2=55;
edge3=77;
binpartition=[0 edge1 edge2 edge3];
agegenehist=[response_array_init(:,3)  response_array_init(:,2)  ]
hist( agegenehist.agegenehist()  );
%}