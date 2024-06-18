function saveFig(figname)
    savefname = char(figname);
    print(savefname,'-dpng');%save as .png
    %print(savefname,'-depsc','-tiff'); %save as .eps level 3 with color, and tiff preview
    %matlabfrag(savefname);
    if ~exist('matlabfrag.m')
        warning('saveFig2:noMatlabFrag','File matlabfrag.m does not exist or is not in MATLAB''s search path.');
        print(savefname,'-depsc','-tiff'); %save as .eps level 3 with color, and tiff preview
    else
        matlabfrag(savefname);
    end

    try
        saveas(gcf,savefname,'fig');% save as .fig
    catch ME1
        idSegLast = regexp(ME1.identifier,'(?<=:)\w+$','match');
        if strcmp(idSegLast,'errorClosingFile') 
            %warning(ME1.message)
            try
                hgsave(gcf,savefname,'-v7.3');
                message=sprintf("%s.fig saved with -v7.3 flag.",figname);
                %warning(figname+".fig saved with -v7.3 flag.")
                warning(char(message))
            catch ME2
                warning(ME2.message)
            end
        end
    end
end