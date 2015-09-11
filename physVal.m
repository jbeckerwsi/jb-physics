classdef physVal < handle
    properties
        value=NaN;
        unit='';
        stdev=NaN;
        UID;
        corelation; %not implemented yet
    end

    
    
    methods
        %constructor
        %val : value
        %stdev: st.deviation
        function obj=physVal(val,stdev,unit)
            obj.value=val;
            if nargin < 3
                unit='';
            end
            if nargin < 2
                stdev=0;
            end
            
            obj.UID=java.rmi.server.UID(); %for corelation-tracking; not implemented yet
            obj.unit=unit;
            obj.stdev=stdev;
            
        end
              
        function r=max(obj1)
           m=-Inf;
           idx=0;
           for i=1:numel(obj1)
               if obj1(i).value>m
                   m=obj1(i).value;
                   idx=i;
               end
           end
           r=obj1(idx);
        end
           
        
        function r = uminus(obj1)
            r=copy(obj1);
            r.value=-r.value;
        end
            
        function r = minus(obj1,obj2)
            if numel(obj1)==numel(obj2)
                for i=1:numel(obj1)
                    r(i)=obj1(i) +(-obj2(i));
                end
            elseif numel(obj2) == 1 && numel(obj1) > 1
                for i=1:numel(obj1)
                    r(i)=obj1(i) + (-obj2);
                end
            end
        end
            
        function r = plus(obj1,obj2)
            if numel(obj1) > 1 && numel(obj2) ==1
                for i=1:numel(obj1)
                    r(i)=obj1(i)+obj2;
                end
                return
            end
            if numel(obj2) > 1 && numel(obj1) ==1
                for i=1:numel(obj2)
                    r(i)=obj1+obj2(i);
                end
                return
            end
            if numel(obj1)==numel(obj2) && numel(obj1) > 1
                for i=1:numel(obj2)
                    r(i)=obj1(i)+obj2(i);
                end
                return
            end
            if isa(obj1,'physVal') && isa(obj2,'physVal')
                if obj1 == obj2 % when calculating stdev for a*b with b=a correlation-coefficient becomes significant ;)
                    r=obj1*2;
                    return
                end
            end
            
            if isa(obj1,'physVal')
                a=obj1.value;
                aErr=obj1.stdev;
                aUnit=obj1.unit;
            else
                a=obj1;
                aErr=0;
                aUnit='';
            end
            
            if isa(obj2,'physVal')
                b=obj2.value;
                bErr=obj2.stdev;
                bUnit=obj2.unit;
            else
                b=obj2;
                bErr=0;
                bUnit='';
            end
            
            
            
            value=a+b;
            stdev=sqrt(aErr^2 + bErr^2);

            if ~strcmp(char(sym(aUnit)),char(sym(bUnit))) && isa(obj1,'physVal') && isa(obj2,'physVal')
                error('addition: units %s and %s dont match!',aUnit,bUnit);
            end
            r=physVal(value,stdev,aUnit);
        end
        
        function r=times(obj1,obj2)
            for i=1:numel(obj1)
                r(i)=obj1(i).mtimes(obj2);
            end
        end
        
        function r=rdivide(obj1,obj2)
            for i=1:numel(obj1)
                r(i)=obj1(i)/obj2;
            end
        end
        
        function r = mtimes(obj1,obj2) %tested 
            if isa(obj1,'physVal') && isa(obj2,'physVal')
                if obj1 == obj2 % when calculating stdev for a*b with b=a correlation-coefficient becomes significant ;)
                    r=obj1^2;
                    return
                end
            end
            
            if isa(obj1,'physVal')
                a=obj1.value;
                aErr=obj1.stdev;
                aUnit=obj1.unit;
            else
                a=obj1;
                aErr=0;
                aUnit='';
            end
            
            if isa(obj2,'physVal')
                b=obj2.value;
                bErr=obj2.stdev;
                bUnit=obj2.unit;
            else
                b=obj2;
                bErr=0;
                bUnit='';
            end
            
            value=a.*b;
            if a==0
                stdev=b.*aErr;
            elseif b==0
                stdev=a.*bErr;
            else
                stdev=value.*sqrt((aErr./a).^2 + (bErr./b).^2);
            end
            
            if strcmp(aUnit,bUnit) && ~strcmp(aUnit,'')
                unit=strcat(aUnit,'^2');
            elseif strcmp(aUnit,'') || strcmp(bUnit,'')
                unit=strcat(aUnit,bUnit);
            else
                unit=strcat(aUnit,'*',bUnit);
            end
            r=physVal(value,stdev,unit);
        end
  
        function r = mrdivide(obj1,obj2)
            if obj1 == obj2 % when calculating stdev for a*b with b=a correlation-coefficient becomes significant ;)
                r=physVal(1,0,'');
                return
            end

            if isa(obj1,'physVal')
                a=obj1.value;
                aErr=obj1.stdev;
                aUnit=obj1.unit;
            else
                a=obj1;
                aErr=0;
                aUnit='1';
            end
            
            if isa(obj2,'physVal')
                b=obj2.value;
                bErr=obj2.stdev;
                bUnit=obj2.unit;
            else
                b=obj2;
                bErr=0;
                bUnit='1';
            end
            
            value=a./b;
            stdev=value.*sqrt((aErr./a).^2 + (bErr./b).^2);
            
            if strcmp(aUnit,bUnit)
                unit='';
            else
                if strcmp(bUnit,'')
                    unit=aUnit;
                else
                    unit=strcat(aUnit,'*((',bUnit,')^-1)');
                end
            end
            r=physVal(value,stdev,unit);
        end
        
        
        function r = mpower(obj,pow)
            if pow == 1
                r=copy(obj);
                return
            end
            
            value=obj.value^pow;
            stdev=value * pow  * (obj.stdev/obj.value);
            unit =sprintf('(%s)^%d',obj.unit,pow);
            r=physVal(value,stdev,unit);
        end
        
       
        function r=copy(obj)
                r=physVal(obj.value,obj.stdev,obj.unit);
        end
            
        function disp(obj)
            if numel(obj) > 1
                for i=1:numel(obj)
                    obj(i).disp();
                end
                return;
            end            
            fprintf('%s\n',obj.str);
        end
        
        function r=abs(obj)
            if numel(obj) > 1
                for i=1:numel(obj)
                    r(i)=abs(obj(i));
                end
                return;
            end
            
            r=physVal(abs(obj.value),obj.stdev,obj.unit);
        end
    
        function s=str(obj,prefix)
            if numel(obj) > 1
                for i=1:numel(obj)
                    s{i}=obj(i).str(prefix);
                end
                return;
            end
           
            
            if nargin < 2
                [factor, prefix]=obj.SiPrefix(obj.value);
            else
                factor=1;
            end
            

            if numel(obj.unit) > 0
                u=char(obj.SIsimplify(char(sym(obj.unit))));
            else
                u='';
            end
            u=strcat(prefix,u);
            s=sprintf('(%.2f %s %.2f)%s',obj.value*factor,char(177),obj.stdev*factor,u);
        end
        

        
        function v=val(obj)
            if numel(obj) == 1
                v=obj.value;
            else
                v=NaN(1,numel(obj));
                for i=1:numel(obj)
                    v(i)=obj(i).value;
                end
            end
        end
        
        
        function c=cell(obj)
            if numel(obj) > 1
                c=cell(1,numel(obj));
                for i=1:numel(obj)
                    c{i}=obj(i);
                end
            else
                c{1}=obj;
            end
        end
        
        function s=sum(obj,~,~)
            s=obj(1);
            if numel(obj) > 1
                for i=2:numel(obj)
                    s=s+obj(i);
                end
            end
        end
    end
    
    
    
    methods (Static=true)
       %physical constants
       %source: http://physics.nist.gov/cgi-bin/cuu/Value?h
       
       function h=h()
            h=physVal(6.62606957E-34,2.9E-41,'J*s'); %Js
       end
       
       function kb=kb()
           %http://physics.nist.gov/cgi-bin/cuu/Value?k|search_for=boltzmann
           kb=physVal(1.3806488E-23, 0.0000013E-23,'J*K^-1');
       end
       
       function m_e=m_e()
           %http://physics.nist.gov/cgi-bin/cuu/Value?me|search_for=electron+mass
           m_e=physVal(9.10938291E-31,0.00000040E-31,'kg');
       end
           
       
       function hbar=hbar()
           hbar=physVal.h / (2*pi);
       end
       
       function e=e()
            e=physVal(1.602176565E-19,3.5E-27,'C'); % C
       end
            
       function e=eps(material)
           if nargin < 1
               material=0;
           end
           switch material
               case 0
                   e=physVal(1/(4*pi*10^-7*299792458^2),0,'F*m^-1'); %nist.gov: e_0 = 1/mu_0*c^2; mu_0=4pi*10^-7
               case 'InAs'
                   e=1;
               case 'HEXtransistor' %O. Wunicke
                   e=2.25;
               otherwise
                   e=1;
           end   
       end
       
       function c=c()
           c=physVal(299792458,0,'m*s^-1');
       end
       
       %%EOF constants
       
       
       function r=acosh(obj)
           if ~strcmp(obj.unit,'')
               warning('acosh: unit not unity');
           end
           %warning('acosh in physVal: error propagation not properly implemented yet');
           r=physVal(acosh(obj.value),0,''); 
       end   
       
       function v=NaN(unit)
           if nargin == 0
               unit='';
           end
           v=physVal(NaN,0,unit);
       end
       
       function r=double(obj)
           if numel(obj) == 1
            r=obj.value;
           else
               for i=1:numel(obj)
                   r(i)=obj.value;
               end
           end
       end
       
       %%% this function expands the unit to SI base units
       function out=SIexpand(input)
           units={'Hz',   'C',   'F',                  'J',          'ohm',             'T',           'V'};
           base ={'s^-1', 'A*s', 's^4*A^2*m^-2*kg^-1', 'kg*m^2*s^-1','kg*m^2*s^-3*A^-2','kg*s^-2*A^-1','kg*m^2*s^-3*A^-1'};
           out=subs(input,units,base);
       end
       
       function out=SIcollapse(input)
           units={'C',   'F',      'F',                  'J',          'ohm',             'T',           'V',               'm^2*V*s',         'Hz'};
           base ={'A*s', 'C*V^-1', 's^4*A^2*m^-2*kg^-1', 'kg*m^2*s^-1','kg*m^2*s^-3*A^-2','kg*s^-2*A^-1','kg*m^2*s^-3*A^-1','kg*m^4*A^-1*s^-2','s^-1'};
           out=subs(input,flip(base),flip(units));
       end  
       
       function out=SIsimplify(input)
           out=physVal.SIcollapse(physVal.SIcollapse(physVal.SIexpand(input)));
       end
       
               % by Jan Simon, Taken from:
        % http://de.mathworks.com/matlabcentral/answers/892-engineering-notation-printed-into-files
        
        function [factor, Str] = SiPrefix(x,greek) % greek: display greek-letters (latex) instead of unicode
            if nargin < 2
                greek=0;
            end
            Exponent = 3 * floor(log10(x) / 3);
            y = x / (10 ^ Exponent);
            ExpValue = [9, 6, 3, 0, -3, -6, -9, -12, -15, -18];
            if greek
                ExpName = {'G', 'M', 'k', '', 'm', '\mu', 'n','p','f','a'};
            else
                ExpName = {'G', 'M', 'k', '', 'm', 'u', 'n','p','f','a'};
            end
            ExpIndex = (Exponent == ExpValue);
            if any(ExpIndex)  % Found in the list:
                Str = ExpName{ExpIndex};
                factor=10^(-Exponent);
            else  % Fallback: Show the numeric value
                % EDITED: Walter refined '%d' to '%+04d':
                Str = '';
                factor=1;
            end
        end
    end
end
    
    