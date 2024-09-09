
import matplotlib.pyplot as plt


def plot_tool_set_axis_properties(cfi, **kwargs):
    """Sets the properties of the axis for a plot.
    
    Args:
        cfi: The ComFiT instance
        kwargs: keyword arguments for the axis properties

    Returns:
        The axis object with the properties set.
    """

    ##### AXIS LIMITS #####
    # xlim is specified as a list
    xlim = [cfi.xmin/cfi.a0, (cfi.xmax-cfi.dx)/cfi.a0]
    if 'xmin' in kwargs:
        xlim[0] = kwargs['xmin'] / cfi.a0
    if 'xmax' in kwargs:
        xlim[1] = kwargs['xmax'] / cfi.a0
    if 'xlim' in kwargs:
        xlim = np.array(kwargs['xlim']) / cfi.a0

    # ylim is specified as a list if dim>1 else as None
    ylim = [cfi.ymin/cfi.a0, (cfi.ymax-cfi.dy)/cfi.a0] if cfi.dim > 1 else None
    if 'ymin' in kwargs:
        ylim[0] = kwargs['ymin'] / cfi.a0 if cfi.dim > 1 else kwargs['ymin']
    if 'ymax' in kwargs:
        ylim[1] = kwargs['ymax'] / cfi.a0 if cfi.dim > 1 else kwargs['ymax']
    if 'ylim' in kwargs:
        ylim = np.array(kwargs['ylim'])/cfi.a0 if cfi.dim > 1 else kwargs['ylim']

    # zlim is specified as a list if dim>2 else as None
    zlim = [cfi.zmin/cfi.a0, (cfi.zmax-cfi.dz)/cfi.a0] if cfi.dim > 2 else None
    if 'zmin' in kwargs:
            zlim[0] = kwargs['zmin'] / cfi.a0 if cfi.dim > 2 else kwargs['zmin']
    if 'zmax' in kwargs:
        zlim[1] = kwargs['zmax'] / cfi.a0 if cfi.dim > 2 else kwargs['zmax']
    if 'zlim' in kwargs:
        zlim = np.array(kwargs['zlim'])/cfi.a0 if cfi.dim > 2 else kwargs['zlim']


    ##### GRID AND TITLE #####
    grid = kwargs.get('grid', True)
    axis_equal = kwargs.get('axis_equal', True)
    title = kwargs.get('title', None)
    suptitle = kwargs.get('suptitle', None)

    ##### SIZE #####
    size = kwargs.get('size', None)

    ##### TICKS #####
    xticks = kwargs.get('xticks', None)  
    xticklabels = kwargs.get('xticklabels', None)

    yticks = kwargs.get('yticks', None)
    yticklabels = kwargs.get('yticklabels', None)

    zticks = kwargs.get('zticks', None)
    zticklabels = kwargs.get('zticklabels', None)

    ##### LABELS #####
    xlabel = kwargs.get('xlabel', 'x/a₀')
    ylabel = kwargs.get('ylabel', 'y/a₀' if cfi.dim > 1 else None)
    zlabel = kwargs.get('zlabel', 'z/a₀' if cfi.dim > 2 else None)

    ##### PLOT NATURE #####
    plot_is_3D = kwargs.get('plot_is_3D', False)

    ##### PLOT LIBRARY #####
    plot_lib = kwargs.get('plot_lib', cfi.plot_lib)
    

    
    ##### AXES #####
    ax = kwargs.get('ax', plt.gca())

    ##### SIZE #####
    if size is not None:
        print("\033[91mWarning: The size keyword is not valid for matplotlib plots.\033[0m")

    ##### TICKS #####
    if xticks is not None:
        ax.set_xticks(xticks)
    if xticklabels is not None:
        ax.set_xticklabels(xticklabels)
    
    if yticks is not None:
        ax.set_yticks(yticks)
    if yticklabels is not None:
        ax.set_yticklabels(yticklabels)

    if zticks is not None:
        ax.set_zticks(zticks)
    if zticklabels is not None:
        ax.set_zticklabels(zticklabels)

    ##### TITLE #####
    if title is not None:
        ax.set_title(title)

    ##### SUPTITLE #####
    if suptitle is not None:
        ax.get_figure().suptitle(suptitle)

    ##### AXIS LABELS #####
    ax.set_xlabel(xlabel)
    
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    
    if zlabel is not None:
        ax.set_zlabel(zlabel)

    ##### AXIS LIMITS #####
    if isinstance(ax, mpl_toolkits.mplot3d.Axes3D):
        ax.set_xlim3d(xlim[0], xlim[1])

        if ylim is not None:
            ax.set_ylim3d(ylim[0], ylim[1])

        if zlim is not None:
            ax.set_zlim3d(zlim[0], zlim[1])

    else:
        ax.set_xlim(xlim[0], xlim[1])

        if ylim is not None:
            ax.set_ylim(ylim[0], ylim[1])

        if zlim is not None:
            ax.set_zlim(zlim[0], zlim[1])

    ##### GRID #####
    ax.grid(grid)

    ##### AXIS ASPECT RATIO #####
    if axis_equal:
        ax.set_aspect('equal')


def plot_field(cfi, field, **kwargs):
        """Plots the given (real) field.
        
    Args:
        field (array-like): The field to be plotted.
        **kwargs: Keyword arguments for the plot.
            See https://comfitlib.com/ClassBaseSystem/ 
            for a full list of keyword arguments.
    
    Returns:
        The axes containing the plot (matplotlib.axes.Axes).

    """

    if field.dtype == bool:
        field = field.astype(float)

    if field.ndim == cfi.dim+1:
        print("\033[91mWarning - in plot_field: the provided field seems to be an order parameter containing several fields. Only the zeroth component will be plotted.\033[0m")
        field = field[0]

    # Check if the vector field is complex
    if np.iscomplexobj(field):
        print("\033[91mWarning - in plot_field: the provided field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
        print('Max imaginary part: ', np.max(np.imag(field)))
        field = np.real(field)

    if plot_lib=='matplotlib':
        # Check if an axis object is provided
        fig = kwargs.get('fig', plt.gcf())
        ax = kwargs.get('ax', None)

    elif plot_lib=='plotly':
        fig = kwargs.get('fig', go.Figure())

    # Kewyord arguments
    colorbar = kwargs.get('colorbar', True)
    axis_equal = kwargs.get('axis_equal', None)

    # Extend the field if not a complete array is given
    field = cfi.plot_tool_extend_field(field)

    kwargs['plot_is_3D'] = False

    ###############################################################
    ###################### DIMENSION: 1 ###########################
    ###############################################################

    if cfi.dim == 1:

        # Keyword arguments particular to the 1D case
        kwargs['grid'] = kwargs.get('grid', True)

        if axis_equal is None:
            kwargs['axis_equal'] = False


        fig.clf()
        ax = fig.add_subplot(111)

        ax.plot(cfi.x/cfi.a0, field)

    ###############################################################
    ###################### DIMENSION: 2 ###########################
    ###############################################################

    if cfi.dim == 2:
        
        # Keyword arguments particular to the 2D case
        kwargs['grid'] = kwargs.get('grid', False)

        ### 2D field Matplotlib ###

        if ax == None:
            fig.clf()
            ax = fig.add_subplot(111)
    
        # Set the colormap
        colormap = kwargs.get('colormap', 'viridis')

        if colormap == 'bluewhitered':
            colormap = tool_colormap_bluewhitered()

        elif colormap == 'sunburst':
            colormap = tool_colormap_sunburst()
        else:
            colormap = plt.get_cmap(colormap)

        # Value limits symmetric
        vlim_symmetric = kwargs.get('vlim_symmetric', False)

        X, Y = np.meshgrid(cfi.x, cfi.y, indexing='ij')

        pcm = ax.pcolormesh(X / cfi.a0, Y / cfi.a0, field, shading='gouraud', cmap=colormap)

        xlim = [cfi.xmin, cfi.xmax-cfi.dx]
        ylim = [cfi.ymin, cfi.ymax-cfi.dy]

        limits_provided = False
        if 'xlim' in kwargs:
            xlim = kwargs['xlim']
            limits_provided = True
        else:
            if 'xmin' in kwargs:
                xlim[0] = kwargs['xmin']
                limits_provided = True
            
            if 'xmax' in kwargs:
                xlim[1] = kwargs['xmax']
                limits_provided = True

        if 'ylim' in kwargs:
            ylim = kwargs['ylim']
            limits_provided = True
        else:
            if 'ymin' in kwargs:
                ylim[0] = kwargs['ymin']
                limits_provided = True
                
            if 'ymax' in kwargs:
                ylim[1] = kwargs['ymax']
                limits_provided = True

        # If explicit limits are provided, use them to change the vlim ranges
        if limits_provided:
            region_to_plot = np.zeros(cfi.dims).astype(bool)
            region_to_plot[(xlim[0] <= X)*(X <= xlim[1])*(ylim[0] <= Y)*(Y <= ylim[1])] = True
            vlim = [np.min(field[region_to_plot]), np.max(field[region_to_plot])]

        else:
            vlim = [np.min(field), np.max(field)]
        
        # Set the value limits
        if 'vlim' in kwargs:
            vlim = kwargs['vlim']
        else:
            if 'vmin' in kwargs:
                vlim[0] = kwargs['vmin']
            if 'vmax' in kwargs:
                vlim[1] = kwargs['vmax']

        if vlim[1] - vlim[0] < 1e-10:
            vlim = [vlim[0]-0.05, vlim[1]+0.05]

        pcm.set_clim(vmin=vlim[0], vmax=vlim[1])

        if 'vlim_symmetric' in kwargs:
            vlim_symmetric = kwargs['vlim_symmetric']
            if vlim_symmetric:
                cmax = abs(field).max()
                cmin = -cmax
                pcm.set_clim(vmin=cmin, vmax=cmax)

        colorbar = kwargs.get('colorbar', True)

        if colorbar:
            cbar = plt.colorbar(pcm, ax=ax)


    ###############################################################
    ###################### DIMENSION: 3 ###########################
    ###############################################################

    elif cfi.dim == 3:

        kwargs['plot_is_3D'] = True

        # Keyword arguments particular to the 3D case

        field_min = np.min(field)
        field_max = np.max(field)

        number_of_layers = kwargs.get('number_of_layers', 1)
        alpha = kwargs.get('alpha', 0.5)
        
        if 'vlim' in kwargs:
            vlim = kwargs['vlim']
            vmin = vlim[0]
            vmax = vlim[1]
        else:
            vmin = field_min
            vmax = field_max

        if 'layer_values' in kwargs:
            layer_values = np.concatenate([[-np.inf], kwargs['layer_values'], [np.inf]])
        else: 
            layer_values = np.linspace(vmin, vmax, number_of_layers + 2)

        if 'colormap' in kwargs:
            colormap = kwargs['colormap']
            if colormap == 'bluewhitered':
                colormap = tool_colormap_bluewhitered()

            elif colormap == 'sunburst':
                colormap = tool_colormap_sunburst()

            else:
                colormap = plt.get_cmap(colormap)
        else: 
            colormap = plt.get_cmap('viridis')
        

        if 'colorbar' in kwargs:
            colorbar = kwargs['colorbar']
        else:
            colorbar = True

        #Plotting the layers

        if ax == None:
            plt.clf()
            ax = plt.gcf().add_subplot(111, projection='3d')


        if field_min < layer_values[1] < field_max:
            cfi.plot_tool_surface_matplotlib(field=field, 
                                    value=layer_values[1], 
                                    color=colormap((layer_values[1]-vmin) / (vmax-vmin)), 
                                    alpha=alpha,
                                    ax=ax)

        for layer_value in layer_values[2:-1]:
            if field_min < layer_value < field_max:
                cfi.plot_tool_surface_matplotlib(field=field, 
                                        value=layer_value, 
                                        color=colormap((layer_value-vmin) / (vmax-vmin)), 
                                        alpha=alpha,
                                        ax=ax)

        if colorbar:
            sm = plt.cm.ScalarMappable(cmap=colormap)
            sm.set_clim(vmin, vmax)
            plt.colorbar(sm, ax=ax, pad=0.2)

    kwargs['ax'] = ax
    ax = plot_tool_set_axis_properties(**kwargs)
    return fig, ax
