from posydon.visualization.VHdiagram import DisplayMode
from posydon.visualization.VH_diagram.Presenter import Presenter, PresenterMode, get_max_distance,CaseInfos,equal_with_epsilon,get_event_state_filename,file_exist,get_star_state_filename
from posydon.visualization.VH_diagram.GraphVisualizer import columnTYPE,StateInfos,ConnectedItem, GraphVisualizerCase

try:
    from PyQt5.QtWidgets import QApplication
    from PyQt5.QtCore import QTimer
except ImportError:
    raise ImportError('PyQt5 is not installed. Please run `pip install .[vis]` in the POSYDON base directory')

import matplotlib.pyplot as plt
from IPython.display import Image, display
from enum import Enum, auto
import os

class VHdiagramm_m:
    """Handle a multiple or hierarchical VH diagram plot."""

    def __init__(
        self,
        filename,
        path="./",
        index=[0],
        *,
        frequency=None,
        hierarchy=False,
        presentMode=PresenterMode.DIAGRAM,
        displayMode=DisplayMode.WINDOW,
        figsize=(10, 8)
    ):
        """Initialize a VHdiagram_m instance."""
        self._app = (
            QApplication.instance()
        )
        if not self._app:
            self._app = QApplication([])
    
        #PRESENTING HIERARCHICALLY
        if hierarchy :
            self._presenter = Presenter_h(filename=filename, path=path, frequency=frequency)
            self._presenter.present_h(index, presentMode)
        #PRESENTING MULTIPLES SIDE BY SIDE
        else:
            self._presenter = Presenter_m(filename=filename, path=path, frequency=frequency)
            self._presenter.present_m(index, presentMode)

        if displayMode == DisplayMode.INLINE_B:
            #void	singleShot(int msec, Functor functor)
            #This static function calls functor after a given time interval.
            QTimer.singleShot(0, lambda: self._display_inline_b())
        elif displayMode == DisplayMode.INLINE_S:
            QTimer.singleShot(0, lambda: self._display_inline_s(figsize))

        self._app.exec_()
    
    #Return indexes witch have been sorted hierarchycally for tree view
    def get_sorted_index(self):
        if type(self._presenter) == Presenter_h :
            return [dataf.index[0] for dataf in self._presenter._columns_data]
        else :#return empty if Presenter_h has no been called
            return []

    def _display_inline_b(self):
        filepath = self._presenter.screen()
        self._presenter.close()
        display(Image(filename=filepath))

    def _display_inline_s(self, figsize):
        filepath = self._presenter.screen()
        self._presenter.close()
        image = plt.imread(filepath)

        plt.figure(figsize=figsize)
        #Display data as an image, i.e., on a 2D regular raster.
        plt.imshow(image, aspect="auto")
        plt.axis("off")
    
class Presenter_m(Presenter):
    """Handle multiple Diagram view"""
    def __init__(self,filename, path="./", frequency=None):
        super(Presenter_m, self).__init__(filename, path)
        self._frequency = frequency
        
    def present_m(self, indexes, mode=PresenterMode.DETAILED):
        """Preset the binary."""
        self._present_mode = mode
        self._visualizer.options().set_showed_mode(self._present_mode)
        self._visualizer().reset()
        
        if self._frequency== None :
            self._frequency = {i: 1/len(indexes) for i in indexes}
    
        for i in indexes :
            self._prepare_corresponding_data(i)
            self._update_visualisation_m()#redefined

        self._main_window.show()

    def _update_visualisation_m(self):
        """Update the display according to current_index and present_mode."""
        if self._current_index is None:
            return

        if self._present_mode == PresenterMode.DETAILED:
            self._prepare_basic_columns_m()
            self._detailed_presentation(self._current_data)
        elif self._present_mode == PresenterMode.REDUCED:
            self._prepare_basic_columns_m()
            self._reduced_presentation(self._current_data)
        elif self._present_mode == PresenterMode.SIMPLIFIED:
            self._prepare_basic_columns_m()
            self._simplified_presentation(self._current_data)
        elif self._present_mode == PresenterMode.DIAGRAM:
            self._prepare_diagram_columns_m()
            self._digram_presentation(self._current_data)
    
    def _prepare_basic_columns_m(self):#THE REST USES 3 COLUMNS COL_SPAN=1
        """Add columns for Detailled/Reduced/Simplified View."""
        self._time_id = self._visualizer().add_column(columnTYPE.TIMELINE)
        self._visualizer().get_column(self._time_id).set_title(f"{self._current_index} : TIME")

        self._S1_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer().get_column(self._S1_id).set_title("S1")

        self._event_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer().get_column(self._event_id).set_title("EVENT/STATE")

        self._S2_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer().get_column(self._S2_id).set_title("S2")
    
    def _prepare_diagram_columns_m(self):
        """Add columns for multiple Diagram View."""
        self._time_id = self._visualizer().add_column(columnTYPE.TIMELINE)
        self._state_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer().get_column(self._state_id).set_title(f"{self._current_index} : {self._frequency[self._current_index] : .2f} %")

class Presenter_h(Presenter):
    """Handles drawing binaries as hierarchical trees"""
    def __init__(self,filename, path="./", frequency=None):
        super(Presenter_h, self).__init__(filename, path)
        self._frequency = frequency
        self._column_width = 400 #fixed but could be dynamic
        
    def present_h(self, indexes, mode=None):
        """Preset the binary."""
        self._present_mode = PresenterMode.DIAGRAM #only available mode
        self._visualizer.options().set_showed_mode(self._present_mode) #keep
        self._visualizer().reset()
        
        if self._frequency == None :#default uniform percentage
            self._frequency = {i: 100/len(indexes) for i in indexes}
        
        self._columns_indexes = []#list of corresponding binaries indexes
        self._columns_data = []#list of corresponding dataframes
        self._infos = []#list of corresponding CaseInfos/StateInfos (with pictures filenames to compare)
        self._edges = []#List of edges in final tree view
    
        for i in indexes :
            #For every index prepare dataframe, add column and record it
            self._prepare_corresponding_data(i)
            self._columns_data.append(self._current_data.copy())
            self._update_visualisation_h()
            
        self._diagram_presentation_h()#sorting removing duplicate setting edges etc
        self._main_window.show()#Graphical présentation
        
    def _sort_h(self):
        """Hierarchical sorting of binaries the way python natureally sort tuples.
        Tuples of représentation filenames here"""
        self._infos.sort(key=lambda info_col : self._tuple_h(info_col))
        #Reajusting column_id witch may be different from when created in _prepare_diagram_line()
        new_col_data = []
        
        for i in range(len(self._infos)):
            new_col_data.append(
                self._columns_data[self._infos[i][0].column_id]
            ) 
            for info in self._infos[i]:                   
                info.column_id = i#don't forget the column_id
        
        self._columns_data = new_col_data

    def _tuple_h(self, info_col):
        """Return a list of filenames of the picture (or centered_text) from the list of 
        StateInfos/CaseInfos for a given column"""
        return [(os.path.split(info.S1_filename)[-1] +
                os.path.split(info.S2_filename)[-1] if (
                    info.event_filename == None
                )else os.path.split(info.event_filename)[-1]) if (
            type(info)==StateInfos
        )else info.centered_text
                for info in info_col]
    
    def _compareInfos(self, info0, info1):
        """Conpare two Infos objects (can be CaseInfos or StateInfos)"""
        if type(info0) == type(info1):
            if type(info0) == StateInfos :
                if info0.event_filename == None == info1.event_filename :
                    return (
                        os.path.split(info0.S1_filename)[-1]==os.path.split(info1.S1_filename)[-1]
                    ) and (
                        os.path.split(info0.S2_filename)[-1]==os.path.split(info1.S2_filename)[-1]
                    )
                else :#One of event_filenames may still be None
                    return info0.event_filename == info1.event_filename
            else :
                return info0.centered_text == info1.centered_text
        else : return False
    
    def _remove_duplicate(self,start_col,
                          end_col,
                          line,
                          parent_col=None):
        """Remove consecutive horizontally duplicated state to get the tree view
        (recursively)"""
        if start_col > end_col :
            return
        if start_col == end_col :
            for info in self._infos[start_col][line : ] :
                #edges are now vertical except maybe the state at line
                info.connected = True
                j = info.column_id
                #binary's relative percentage
                rel = self._frequency[self._columns_data[j].index[0]]
                if type(info)==StateInfos :
                    info.bot_texts = [f"{self._columns_data[start_col].index[0]} : {rel: .30f} %"]
                else :
                    info.bot_left_text = f"{self._columns_data[start_col].index[0]} : {rel: .30f} %"
                    
            if (0 < line < len(self._infos[start_col]) ) and (
                type(self._infos[start_col][line-1]) == CaseInfos
            ) :
                self._infos[start_col][line].connected = False
                self._edges.extend([[line-1,
                             parent_col,
                             line, 
                             start_col] ])
            
            return

        shortest_col = min(self._infos[start_col : end_col+1] ,
                       key=lambda info_col: len(info_col))
        
        #if no column is shorter that the length of the line
        if len(shortest_col) > line:
            info_col0 = self._infos[start_col]
            #Range of columns with same state at current line
            rg_column = [[start_col]]
            
            for j in range(start_col+1, end_col+1):
                info_col1 = self._infos[j]
                if self._compareInfos(info_col0[line], info_col1[line]):
                    #replacing by empty state
                    info_col0[line] = CaseInfos(
                        info_col0[line].column_id,
                        centered_txt = " "
                    )
                    info_col0[line].connected = False
                else:
                    #recording a new range
                    rg_column[-1].append(j-1)
                    rg_column.append([j])
                #righ state become comparaison basis for the loop
                info_col0 = info_col1
            rg_column[-1].append(end_col)
            #an edge (A --> B) is represented by [x1, y1, x2, y2] where A:(x1,y1) and B:(x2,y2)
            self._edges.extend([[line-1,
                                 parent_col,
                                 line, 
                                 c[0] + (c[1]-c[0])//2] for c in rg_column])
            
            self._center_column(rg_column, line)
            
            #recursive call with line+1
            for start,end in rg_column :
                self._remove_duplicate(start, end, line+1,
                                      parent_col = start + (end - start)//2)
        
        #if one column is shorter than the line's length separate range in two
        else :#(there may be more than one shorter column and/or in sequence)
            no_col = self._infos.index(shortest_col, start_col, end_col+1)
            
            self._remove_duplicate(start_col, no_col-1, line, parent_col = parent_col)
            self._remove_duplicate(no_col+1, end_col, line, parent_col = parent_col)
            
        
    def _center_column(self, range_column, line):
        """move states in a range of columns from far right to center
        also setting the sum of same states's frequencies appearance to be printed"""
        for k,j in range_column:
            if j - k >= 1 :
                info = self._infos[j][line]
                self._infos[j][line] = CaseInfos(
                        info.column_id,
                        centered_txt = " "
                    )
                self._infos[j][line].connected = False
                
                info.column_id = self._infos[k + (j-k)//2][line].column_id
                relative_freq = sum(
                    [self._frequency[self._columns_data[i].index[0]]
                    for i in range(k,j+1)]
                )
                if type(info)==StateInfos :
                    info.bot_texts = [f"{self._columns_data[k].index[0]} to {self._columns_data[j].index[0]} : {relative_freq: .25f} %"]
                else :
                    info.bot_left_text = f"{self._columns_data[k].index[0]} to {self._columns_data[j].index[0]} : {relative_freq: .25f} %"
                #swapping
                self._infos[k + (j-k)//2][line] = info
            elif j==k :
                info = self._infos[j][line]
                relative_freq = self._frequency[self._columns_data[j].index[0]]
                if type(info)==StateInfos :
                    info.bot_texts = [f"{self._columns_data[k].index[0]} to {self._columns_data[j].index[0]} : {relative_freq: .25f} %"]
                else :
                    info.bot_left_text = f"{self._columns_data[k].index[0]} to {self._columns_data[j].index[0]} : {relative_freq: .25f} %"      
        

    def _update_visualisation_h(self):
        """Hierarchy PresenterMode for the moment is only PresenterMode.DIAGRAM"""
        if self._current_index is None:
            return

        if self._present_mode == PresenterMode.DETAILED:
            self._prepare_basic_columns()#irrelevant
        elif self._present_mode == PresenterMode.REDUCED:
            self._prepare_basic_columns()#irrelevant
            self._reduced_presentation(self._current_data)
        elif self._present_mode == PresenterMode.SIMPLIFIED:
            self._prepare_basic_columns()#irrelevant
            self._simplified_presentation(self._current_data)
        elif self._present_mode == PresenterMode.DIAGRAM:
            self._prepare_diagram_columns_h()
            
    
    def _prepare_diagram_columns_h(self):
        """Add columns for Diagram view set column minimum width and add index to _columns_indexes"""
        self._state_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer()._layout.setColumnMinimumWidth(self._state_id, self._column_width)
        self._columns_indexes.append(self._state_id)
    
    def _diagram_presentation_h(self):
        
        simplified_datas = []#list of columns's simplified datas
        for i in self._columns_indexes:
            simplified_datas.append(
                self._simplify_data(self._columns_data[i])
            )
        
        for i in self._columns_indexes:
            simplified_data = simplified_datas[i]
            self._state_id = i
            max_distance_data = max(simplified_data, key=get_max_distance)
            max_distance = get_max_distance(max_distance_data)

            self._visualizer()._columns[i]._row_index = 16#forget why
            
            #Filling the list of StateInfos of each lines of the column
            self._infos.append([])
            for line_data in simplified_data:
                dline = self._prepare_diagram_line(line_data, max_distance)
                self._infos[-1].append(dline[0])
            
            #We add a CaseInfos according to state_after of last simplified_datas
            if (
                simplified_data[-1]["state"].state_after == "disrupted"
                or simplified_data[-1]["state"].state_after == "merged"
            ):
                aditional_info = CaseInfos(self._state_id)
                aditional_info.border_width = 2
                aditional_info.centered_text = simplified_data[-1][
                    "state"].state_after
                self._infos[-1].append(aditional_info)
                aditional_info.connected = False
                
        self._sort_h()#Sorting columns hierarchically according to pictures filenames
        self._remove_duplicate(0,len(self._infos)-1,0)#Removing horizontal duplicate 
        
        #Setting percentage of each column to be printed (above the column)
        for i in self._columns_indexes:
            self._state_id = i
            binary_index = self._columns_data[i].index[0]
            self._visualizer()._columns[i]._row_index = 16 #forget why
            
            self._visualizer().get_column(i).set_title(f"{binary_index} : {self._frequency[binary_index]: .4f} %")
            for line_info in self._infos[i]:
                self._visualizer().add_line([line_info])
           
        self._paint_obl_edges()
    
    def _prepare_diagram_line(self, data, max_distance):
        
        state_info = StateInfos(self._state_id)
        if self._distance_representation and "separation" in data:
            state_info.distance = self._get_distance_representation(
                data["separation"], max_distance
            )
        else:
            state_info.distance = 1

        state_before_filename = os.path.join(
            self.PATH_TO_DRAWS,
            get_event_state_filename(
                data["S1_state"].state_before,
                data["state"].state_before,
                data["S2_state"].state_before,
                suffix=".png",
            ),
        )
        state_after_filename = os.path.join(
            self.PATH_TO_DRAWS,
            get_event_state_filename(
                data["S1_state"].state_before,
                data["state"].state_before,
                data["S2_state"].state_before,
                suffix=".png",
            ),
        )
        if (data["state"].state_after != "detached"
                and file_exist(state_after_filename)):
            state_info.event_filename = state_after_filename
        elif data["state"].state_before != "detached" and file_exist(
            state_before_filename
        ):
            state_info.event_filename = state_before_filename
        else:
            if data["S1_state"].state_after is not None:
                state_info.S1_filename = os.path.join(
                    self.PATH_TO_DRAWS,
                    get_star_state_filename(
                        data["S1_state"].state_after, suffix=".png"
                    ),
                )
            else:
                state_info.S1_filename = os.path.join(
                    self.PATH_TO_DRAWS,
                    get_star_state_filename(
                        data["S1_state"].state_before, suffix=".png"
                    ),
                )

            if data["S2_state"].state_after is not None:
                state_info.S2_filename = os.path.join(
                    self.PATH_TO_DRAWS,
                    get_star_state_filename(
                        data["S2_state"].state_after, suffix=".png"
                    ),
                )
            else:
                state_info.S2_filename = os.path.join(
                    self.PATH_TO_DRAWS,
                    get_star_state_filename(
                        data["S2_state"].state_before, suffix=".png"
                    ),
                )
        state_info.connected = False
        #Traditionnaly returning a list with StateInfos/CaseInfos/PointInfos so we return a list  
        return [state_info]
    
    def _paint_obl_edges(self):
        """Adding oblic edges to be painted as ConnectedItem from edges added in self._edges"""
        for i1,j1,i2,j2 in self._edges[1:] :
            self._visualizer()._columns[j1]._connected_items.append(
                        ConnectedItem(
                            self._visualizer()._columns[j1]._items[i1],
                            self._visualizer()._columns[j2]._items[i2]))