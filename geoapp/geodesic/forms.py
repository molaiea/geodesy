from django import forms
from django.core.validators import MaxValueValidator
from django.forms import MultiWidget, NumberInput

<<<<<<< HEAD


=======
>>>>>>> 4379fda8acb039391c8ec2bfe2ad5282d000b36b
class LatitudeWidget(forms.MultiWidget):
    def __init__(self, attrs=None):
        super().__init__([
            forms.NumberInput(attrs={'max':90, 'min':0}),
            forms.NumberInput(attrs={'max':59, 'min':0}),
            forms.NumberInput(attrs={'max':59, 'min':0}),
            forms.Select(choices=(
            ('north', 'N'),
            ('south', 'S')))])
        
    
    def decompress(self, value):
        if value:
            return value
        return [0, 0, 0, 'N']

class LatitudeField(forms.MultiValueField):
    widget = LatitudeWidget

    def __init__(self, *args, **kwargs):
        fields = (
            forms.IntegerField(),
            forms.IntegerField(),
            forms.IntegerField(),
            forms.ChoiceField(choices=[('north', 'N'),('south', 'S')])
        )

        super().__init__(fields, *args, **kwargs)
    
    def compress(self, data_list):
        if data_list[3] == 'north':
            return data_list[0]+data_list[1]/60+data_list[2]/3600
        return -(data_list[0]+data_list[1]/60+data_list[2]/3600)

class LongitudeWidget(forms.MultiWidget):
    def __init__(self, attrs=None):
        super().__init__([
            forms.NumberInput(attrs={'max':180, 'min':0}),
            forms.NumberInput(attrs={'max':59, 'min':0}),
            forms.NumberInput(attrs={'max':59, 'min':0}),
            forms.Select(choices=(
            ('est', 'E'),
            ('west', 'O')))])
        
    
    def decompress(self, value):
        if value:
            return value
        return [0, 0, 0, 'E']

class LongitudeField(forms.MultiValueField):
    widget = LongitudeWidget

    def __init__(self, *args, **kwargs):
        fields = (
            forms.IntegerField(),
            forms.IntegerField(),
            forms.IntegerField(),
            forms.ChoiceField(choices=[('est', 'E'),('west', 'O')])
        )

        super().__init__(fields, *args, **kwargs)
    
    def compress(self, data_list):
        if data_list[3] == 'est':
            return data_list[0]+data_list[1]/60+data_list[2]/3600
        return -(data_list[0]+data_list[1]/60+data_list[2]/3600)

class parametersWidget(forms.MultiWidget):
    def __init__(self, attrs=None):
        super().__init__([
            forms.NumberInput(),
            forms.NumberInput()])
        
    
    def decompress(self, value):
        if value:
            return value
        return [0, 0]

class parametersField(forms.MultiValueField):
    widget = parametersWidget

    def __init__(self, *args, **kwargs):
        fields = ([
            forms.FloatField(),
            forms.FloatField()])
        

        super().__init__(fields, *args, **kwargs)
    
    def compress(self, data_list):
        if data_list:
            return data_list
        return [None, None]
class finalform(forms.Form):
    ellipsoid = forms.ChoiceField(choices=[('wgs', 'WGS84'),('grs', 'GRS80'), ('clarke', 'Clarke')], required=False)
    params =  parametersField(required=False)
    latitude = LatitudeField()
    longitude = LongitudeField()
    azimut = forms.FloatField()
    distance_geodesique = forms.FloatField()
    