classdef plot_save < handle
    
    properties
        x_width;
        y_width;
        y_width2;
    end
    
    methods
        function obj = plot_save(x,y,y2)
            obj.x_width=x;
            obj.y_width=y;
            obj.y_width2=y2;
        end
        
        function save_eps(obj,str)
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperPosition', [0 0 obj.x_width obj.y_width]);
            str=strcat('graphs\',str);
            strfig=strcat(str,'.fig');
            savefig(strfig);
            streps=strcat(str,'.eps');
            saveas(gcf,streps,'epsc')
        end

        function save_png(obj,str)
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperPosition', [0 0 obj.x_width obj.y_width]);
            str=strcat('graphs\',str);
            strfig=strcat(str,'.fig');
            savefig(strfig);
            strpng=strcat(str,'.png');
            saveas(gcf,strpng,'png')
        end

        function save_png2(obj,str)
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperPosition', [0 0 obj.x_width obj.y_width2]);
            str=strcat('graphs\',str);
            strfig=strcat(str,'.fig');
            savefig(strfig);
            strpng=strcat(str,'.png');
            saveas(gcf,strpng,'png')
        end
        
    end
end