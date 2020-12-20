# ----------------------------------------------------------------------------------------------------------------------
struct CoordsToIndices
    XTopLeft     :: Float64
    YTopLeft     :: Float64
    XBottomRight :: Float64
    YBottomRight :: Float64
    Width        :: Int64
    Height       :: Int64
    function CoordsToIndices(coord1:: Tuple{Number, Number},
                             coord2:: Tuple{Number, Number},
                             width :: Int64,
                             height:: Int64)
        X1, Y1 = coord1
        X2, Y2 = coord2
        if !(X1 < X2 && Y1 > Y2)
            throw("Coordinate 1 is the top left, and Coordinate 2 is the bottom right.")
        end
        if !(width >= 1 && height >= 1)
            throw("The height and width of the window has to be larger than 1 by 1")
        end
        new(X1, Y1, X2, Y2, width, height)
    end
end

function GetImageIndices(obj::CoordsToIndices, x, y)
    HeightInCoords  = obj.YTopLeft - obj.YBottomRight
    WidthInCoords   = obj.XBottomRight - obj.XTopLeft
    OffSetXInCoords = x - obj.XTopLeft
    OffSetYInCoords = obj.YTopLeft - y
    OffSetXRatio    = OffSetXInCoords/WidthInCoords
    OffSetYRatio    = OffSetYInCoords/HeightInCoords
    OffSetXPixel    = ceil(OffSetXRatio*obj.Width)
    OffSetYPixel    = ceil(OffSetYRatio*obj.Height)
    return trunc(Int64, OffSetYPixel), trunc(Int64, OffSetXPixel)
end

# ----------------------------------------------------------------------------------------------------------------------