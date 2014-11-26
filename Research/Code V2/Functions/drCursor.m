function drCursor(C,type,sty)
if ~exist('type','var')
    type = 'Vertical';
end
if ~exist('sty','var')
    sty = 'r-';
end

NC = numel(C);
if(size(type,1) == 1)
    type = repmat(type,NC,1);
end
if(size(sty,1) == 1)
    sty = repmat(sty,NC,1);
end

AX = gca;
hold on;
XLIM = get(AX,'XLim');
YLIM = get(AX,'YLim');

for i=1:NC
    switch lower(type(i,:))
        case 'vertical'
           x0 = C(i); 
           x1 = C(i);
           y0 = YLIM(1);
           y1 = YLIM(2);
        case 'horizontal'
           y0 = C(i); 
           y1 = C(i);
           x0 = XLIM(1);
           x1 = XLIM(2);
        otherwise 
            error('Cursor type must be either ''Vertical'' or ''Horizontal''.');
    end
    plot([x0 x1], [y0 y1], sty(i,:));
end
end